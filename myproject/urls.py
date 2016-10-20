"""myproject URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import include, url
from django.contrib import admin
#from myapp.views import *
from myapp import views
from myapp import urls

from django.conf import settings
from django.conf.urls.static import static
from django.views.generic import RedirectView

urlpatterns = [
    url(r'^/', include('myapp.urls',namespace='myapp')),

    url(r'^admin/', admin.site.urls),
    url(r'^$',views.home),
    url(r'^datasets/', views.datasets),
    url(r'^metrics/', views.metrics),
    url(r'^python_code/', views.python_code),
    url(r'^docs/', views.docs),
    url(r'^frequent_questions/', views.frequent_questions),
    url(r'^email/', views.email),
    url(r'^output/', views.output),

    url(r'^upload/', views.upload),
    url(r'^process_sort/', views.process_sort),

    #url(r'^polls/', include('polls.urls', namespace="polls")),

    #url(r'^home/$', include('myapp.urls',namespace='myapp')),

    #url(r'^list_dir/$', views.list_dir, name="list_dir/"), #Can't link directly as .html file needs data from the views
    #url(r'^simple_chart/$', views.simple_chart, name="simple_chart"),
    #url(r'^/$', simple_chart, name="simple_chart"),
    
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
