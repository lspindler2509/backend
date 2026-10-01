import base64
import datetime
import json
import os
import random
import string
import sys
import time

import requests

from drugstone.settings import DEFAULTS
from drugstone.util.property_calulations import calculate_properties_id_based
from tasks.task_hook import TaskHook
from tasks.util.custom_network import add_edges, filter_proteins, remove_ppi_edges
from tasks.util.read_graph_tool_graph import read_graph_tool_graph

# Base URL for KeyPathwayMiner Web API
KPM_URL = 'https://exbio.wzw.tum.de/keypathwayminer/requests/'


def kpm_task(task_hook: TaskHook):
    """
    Run KeyPathwayMiner using the currently active DrugstOne network.
    Uploads the active PPI graph as a custom network to KeyPathwayMinerWeb via multipart/form-data,
    polls progress asynchronously, and formats results back into DrugstOne's standard network structure.

    :param task_hook: TaskHook instance with seeds, config, and dataset parameters.
    """
    config = task_hook.parameters.get("config", {})
    id_space = config.get("identifier", "symbol")
    is_reviewed = config.get("reviewed", False)

    ppi_dataset = task_hook.parameters.get("ppi_dataset")
    if not ppi_dataset:
        ppi_dataset = {"name": DEFAULTS.get("ppi", "NeDRex"), "licenced": False}

    pdi_dataset = task_hook.parameters.get("pdi_dataset")
    if not pdi_dataset:
        pdi_dataset = {"name": DEFAULTS.get("pdi", "NeDRex"), "licenced": False}

    seeds = list(task_hook.seeds)

    # --- 1. Load the active network
    task_hook.set_progress(0.05, "Loading active network")

    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset.get("licenced") or pdi_dataset.get("licenced"):
        filename += "_licenced"
    if is_reviewed:
        filename += "_reviewed"
    file_path = os.path.join(task_hook.data_directory, filename + ".gt")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Network file not found: {file_path}")

    # Read graph, target='protein' strips drug nodes and keeps only the PPI network
    g, seed_ids, _ = read_graph_tool_graph(file_path, seeds, id_space, sys.maxsize, target="protein")

    # Apply custom edges and/or custom nodes if configured
    custom_edges = task_hook.parameters.get("custom_edges", False)
    no_default_edges = task_hook.parameters.get("exclude_drugstone_ppi_edges", False)
    custom_nodes = task_hook.parameters.get("network_nodes", False)

    if custom_edges:
        if no_default_edges:
            g = remove_ppi_edges(g)
        edges = task_hook.parameters.get("input_network", {}).get("edges", [])
        g = add_edges(g, edges)

    if custom_nodes:
        g, seed_ids, _ = filter_proteins(g, custom_nodes, [], seeds)

    # --- 2. Serialize graph to 2-column TSV edge list
    task_hook.set_progress(0.1, "Serializing network for KeyPathwayMiner")

    node_attr = "internal_id"
    edges_seen = set()
    network_lines = []

    for e in g.edges():
        u = g.vertex_properties[node_attr][e.source()]
        v = g.vertex_properties[node_attr][e.target()]
        if not u or not v or u == v:
            continue
        edge_key = (min(u, v), max(u, v))
        if edge_key not in edges_seen:
            edges_seen.add(edge_key)
            network_lines.append(f"{u}\t{v}\n")

    if not network_lines:
        raise RuntimeError("No edges found in the active network to send to KeyPathwayMiner.")

    network_tsv = "".join(network_lines)

    # --- 3. Build indicator matrix with active seeds
    task_hook.set_progress(0.15, "Preparing indicator matrix")

    indicator_tsv = "".join(f"{seed}\t1\n" for seed in seeds)
    content_b64 = base64.b64encode(indicator_tsv.encode("utf-8")).decode("ascii")

    attached_to_id = "".join(random.choices(string.ascii_uppercase + string.digits, k=32))
    dataset_name = "indicatorMatrix"

    datasets = [
        {
            "name": dataset_name,
            "fileName": "indicator.tsv",
            "attachedToID": attached_to_id,
            "hasHeader": False,
            "valueType": "binary",
            "contentBase64": content_b64,
        }
    ]

    # --- 4. Configure KPM settings (graphID is omitted so KPM uses the custom uploaded graph)
    k_val = int(task_hook.parameters.get("k", 1))
    computed_pathways = int(task_hook.parameters.get("computed_pathways", 1))

    kpm_settings = {
        "parameters": {
            "name": f"Drugstone run on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "algorithm": "GREEDY",
            "strategy": "INES",
            "removeBENs": "true",
            "unmapped_nodes": "Add to negative list",
            "computed_pathways": computed_pathways,
            "k_values": {
                "val": k_val,
            },
            "l_values": [
                {
                    "val": 0,
                    "datasetName": dataset_name,
                }
            ],
        },
        "withPerturbation": "false",
        "linkType": "OR",
        "attachedToID": attached_to_id,
    }

    # --- 5. Submit job asynchronously via multipart/form-data
    task_hook.set_progress(0.2, "Submitting analysis to KeyPathwayMiner")

    submit_url = KPM_URL + "submitAsync"
    payload_data = {
        "kpmSettings": json.dumps(kpm_settings),
        "datasets": json.dumps(datasets),
        "networkFileName": "network.tsv",
        "networkHasHeader": "false",
    }
    files = {
        "graphFile": ("network.tsv", network_tsv.encode("utf-8"), "text/tab-separated-values"),
    }

    try:
        response = requests.post(submit_url, data=payload_data, files=files, timeout=60)
        response.raise_for_status()
        submit_json = response.json()
    except Exception as e:
        raise RuntimeError(f"Failed to submit job to KeyPathwayMiner: {e}")

    if not submit_json.get("success"):
        raise RuntimeError(f"Job submission failed. Server response:\n{submit_json}")

    quest_id = submit_json["questID"]

    # --- 6. Poll run status until completed
    task_hook.set_progress(0.25, "Queued in KeyPathwayMiner")
    status_url = KPM_URL + f"runStatus?questID={quest_id}"

    old_progress = -1
    while True:
        try:
            status_resp = requests.get(status_url, timeout=15)
            status_resp.raise_for_status()
            status_json = status_resp.json()
        except Exception:
            time.sleep(1)
            continue

        if not status_json.get("runExists", True):
            raise RuntimeError(f"Job status retrieval failed. Run does not exist:\n{status_json}")

        if status_json.get("failed"):
            raise RuntimeError(f"KeyPathwayMiner run failed: {status_json.get('statusMessage', 'Unknown error')}")

        progress = float(status_json.get("progress", 0.0))
        scaled_progress = 0.25 + (progress * 0.65)
        if scaled_progress != old_progress:
            status_msg = status_json.get("statusMessage", "Running KeyPathwayMiner...")
            task_hook.set_progress(progress=scaled_progress, status=status_msg)
            old_progress = scaled_progress

        if status_json.get("completed") or status_json.get("cancelled"):
            break

        time.sleep(1)

    if status_json.get("cancelled"):
        raise RuntimeError("KeyPathwayMiner run was cancelled.")

    # --- 7. Retrieve and parse results
    task_hook.set_progress(0.9, "Retrieving results from KeyPathwayMiner")
    results_url = KPM_URL + f"results?questID={quest_id}"

    try:
        results_resp = requests.get(results_url, timeout=30)
        results_resp.raise_for_status()
        results_json = results_resp.json()
    except Exception as e:
        raise RuntimeError(f"Failed to retrieve results from KeyPathwayMiner: {e}")

    if not results_json.get("success"):
        raise RuntimeError(f"KeyPathwayMiner completed but was unsuccessful:\n{results_json}")

    returned_nodes = set()
    returned_edges = set()

    if results_json.get("subnetworks"):
        for sub_id, sub_data in results_json["subnetworks"].items():
            for node in sub_data.get("nodes", []):
                node_name = node if isinstance(node, str) else (node.get("name") or node.get("id"))
                if node_name:
                    returned_nodes.add(node_name)
            for edge in sub_data.get("edges", []):
                s = edge.get("source")
                t = edge.get("target")
                if s and t:
                    returned_edges.add((min(s, t), max(s, t)))
    elif results_json.get("union_network"):
        u = results_json["union_network"]
        for node in u.get("nodes", []):
            node_name = node if isinstance(node, str) else (node.get("name") or node.get("id"))
            if node_name:
                returned_nodes.add(node_name)
        for edge in u.get("edges", []):
            s = edge.get("source")
            t = edge.get("target")
            if s and t:
                returned_edges.add((min(s, t), max(s, t)))
    elif results_json.get("resultGraphs"):
        for graph in results_json["resultGraphs"]:
            if graph.get("isUnionSet"):
                continue
            for node in graph.get("nodes", []):
                node_name = node.get("name") or node.get("id")
                if node_name:
                    returned_nodes.add(node_name)
            for edge in graph.get("edges", []):
                s = edge.get("source")
                t = edge.get("target")
                if s and t:
                    returned_edges.add((min(s, t), max(s, t)))

    # --- 8. Format results for DrugstOne frontend
    task_hook.set_progress(0.95, "Formatting results")

    accepted_nodes = sorted(list(returned_nodes))
    edges_unique = [{"from": s, "to": t} for s, t in sorted(list(returned_edges))]

    subgraph = {
        "nodes": accepted_nodes,
        "edges": edges_unique,
    }

    seeds_set = set(seeds)
    node_types = {node: "protein" for node in accepted_nodes}
    is_seed = {node: (node in seeds_set) for node in accepted_nodes}
    target_nodes = [node for node in accepted_nodes if node not in seeds_set]

    calculateProperties = config.get("calculate_properties", False)
    properties = calculate_properties_id_based(accepted_nodes, g, edges_unique, calculateProperties)

    result_dict = {
        "network": subgraph,
        "target_nodes": target_nodes,
        "node_attributes": {"node_types": node_types, "is_seed": is_seed},
        "properties": properties,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
    }
    task_hook.set_results(results=result_dict)
