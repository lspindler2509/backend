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

from drugstone.views import (
    map_nodes,
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
)

# cache time is 6 hours
urlpatterns = [
    path("get_datasets/", get_datasets),
    path("map_nodes/", map_nodes),
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
    path("get_license", get_license),
    path("send_bugreport/", send_bugreport),
    path("save_selection", save_selection),
    path("view/", get_view),
    path("view_infos", get_view_infos),
]
