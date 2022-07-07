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

from drugstone.views import map_nodes, tasks_view, result_view, \
    graph_export, query_proteins, TissueView, TissueExpressionView, query_tissue_proteins, TaskView, \
    adjacent_drugs, adjacent_disorders, fetch_edges, create_network, load_network

# cache time is 6 hours
urlpatterns = [
    # path('network/', cache_page(21600)(ProteinViralInteractionView.as_view())),
    path('map_nodes/', map_nodes),
    path('fetch_edges/', fetch_edges),
    path('task/', TaskView.as_view()),
    path('tasks/', tasks_view),
    path('task_result/', result_view),
    path('graph_export/', graph_export),
    path('query_proteins/', query_proteins),
    path('query_tissue_proteins/', query_tissue_proteins),
    path('adjacent_drugs/', adjacent_drugs),
    path('adjacent_disorders/', adjacent_disorders),
    # path('drug_interactions/', ProteinDrugInteractionView.as_view()),
    path('tissue_expression/', TissueExpressionView.as_view()),
    path('tissues/', TissueView.as_view()),
    path('admin/', admin.site.urls),
    path('create_network', create_network),
    path('load_network', load_network)
]
