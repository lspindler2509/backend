"""drugstone URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.urls import re_path

from drugstone.views import (
    FileUploadView,
    download_network,
    get_all_scores_pathway_enrichment,
    get_pathway_sources,
    map_nodes,
    overlay_directed_edges,
    prepare_pruning,
    prune,
    recalculate_statistics,
    rename_selection,
    tasks_view,
    result_view,
    graph_export,
    TissueView,
    TissueExpressionView,
    query_tissue_proteins,
    TaskView,
    adjacent_drugs,
    adjacent_disorders,
    fetch_edges,
    create_network,
    load_network,
    get_license,
    get_datasets,
    get_max_tissue_expression,
    convert_compact_ids,
    get_default_params,
    send_bugreport,
    save_selection,
    get_view,
    get_view_infos,
    calculate_result_for_pathway,
    create_genesets,
    add_edges,
    apply_layout,
    searchProteins,
    update_network,
rename_task
)

# cache time is 6 hours
urlpatterns = [
    path("get_datasets/", get_datasets),
    path("download_network", download_network),
    path("map_nodes/", map_nodes),
    path("apply_layout/", apply_layout),
    path("convert_compact_node_list/", convert_compact_ids),
    path("fetch_edges/", fetch_edges),
    path("task/", TaskView.as_view()),
    path("tasks/", tasks_view),
    path("task_result/", result_view),
    path("graph_export/", graph_export),
    path("query_tissue_proteins/", query_tissue_proteins),
    path("adjacent_drugs/", adjacent_drugs),
    path("adjacent_disorders/", adjacent_disorders),
    path("tissue_expression/", TissueExpressionView.as_view()),
    path("tissue_max_expression/", get_max_tissue_expression),
    path("tissues/", TissueView.as_view()),
    path("admin/", admin.site.urls),
    path("create_network", create_network),
    path("load_network", load_network),
    path("get_default_params", get_default_params),
    path("get_pathway_sources", get_pathway_sources),
    path("get_license", get_license),
    path("send_bugreport/", send_bugreport),
    path("save_selection", save_selection),
    path("rename_selection", rename_selection),
    path("rename_task", rename_task),
    path("view/", get_view),
    path("view_infos", get_view_infos),
    path("calculate_result_for_pathway/", calculate_result_for_pathway),
    path("get_all_scores_pathway_enrichment/", get_all_scores_pathway_enrichment),
    path("create_genesets/", create_genesets),
    path("add_edges/", add_edges),
    path("search_proteins/", searchProteins),
    path("prepare_pruning/", prepare_pruning),
    path("prune/", prune),
    path("update_result_network/", update_network),
    path("recalculate_statistics/", recalculate_statistics),
    path("overlay_directed_edges/", overlay_directed_edges),
    re_path(r'^upload/(?P<filename>[^/]+)$', FileUploadView.as_view())

]
