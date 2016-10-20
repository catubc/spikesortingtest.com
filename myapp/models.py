from __future__ import unicode_literals

from django.db import models

from django.conf import settings
from django.core.files.storage import FileSystemStorage
# from data_importer.importers import CSVImporter

class Document(models.Model):
    #docfile = models.FileField(upload_to='uploads/%Y/%m/%d')
    docfile = models.FileField(upload_to='uploads/')
    
    def __unicode__(self):
        return '%s' % (self.docfile.name)
        
    def delete(self, *args, **kwargs):
        os.remove(os.path.join(settings.MEDIA_ROOT, self.docfile.name))
        super(Document,self).delete(*args,**kwargs)

class Result(models.Model):
    dummy = 1


# Create your models here.
