from django.conf.urls import include, url
from django.conf import settings
from django.conf.urls.static import static
from django.views.generic import RedirectView

from django.contrib import admin
from . import views

#app_name = 'myapp'
urlpatterns = [
    #url('', views.list_dir,name='list_dir'),
    #url('', views.hello, name='hello'),
    #url(r'^$', views.list_dir, name='list_dir'),
    url(r'datasets/', views.datasets, name='datasets'),
    url(r'upload/', views.upload, name='upload'),
    url(r'metrics/', views.metrics, name='metrics'),
    url(r'python_code/', views.python_code, name='python_code'),
    url(r'docs/', views.docs, name='docs'),
    url(r'frequent_questions/', views.frequent_questions, name='frequent_questions'),
    url(r'email/', views.email, name='email'),
    url(r'output/', views.output, name='output'),
    url(r'process_sort/', views.process_sort, name='process_sort'),
    


    #url(r'^home/$', views.home, name='home'),
    #url(r'^home/$', include('myapp.urls',namespace='myapp')),

    #url('^simple_chart$', views.simple_chart,name='simple_chart'),
] 
