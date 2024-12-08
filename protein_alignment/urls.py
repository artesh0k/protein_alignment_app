from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('sub_mat/', views.substitution_matrices_info, name='sub_mat'),
    path('get_alignment_results/', views.get_alignment_results, name='get_alignment_results'),

]
