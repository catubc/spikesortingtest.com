from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader, Context
from django.core.mail.message import EmailMessage
from django.shortcuts import redirect

#Custom script
import pdb
import csv
# -*- coding: utf-8 -*-
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.conf import settings

#from .models import Document
from .forms import DocumentForm, ContactForm
from .models import Document

#Django-importer
#import data_importer
#from data_importer.importers import CSVImporter
#from .models import MyCSVImporterModel

from django.shortcuts import render
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import components, autoload_static

from . import dirk_2
from dirk_2 import *

# From server example (ImportError: No module named 'bokeh.session')
# http://bokeh.pydata.org/en/latest/docs/user_guide/embed.html
from bokeh.plotting import figure, push
from bokeh.embed import autoload_server
# from bokeh.session import Session
#from bokeh.document import Document

# From csv upload: http://stackoverflow.com/questions/14160221/django-read-upload-a-csv
from django.template import RequestContext
from django.shortcuts import render_to_response
# from web.forms import codeUploadForm
# from web.csvTools import CodeCSvModel

def hello(request):
    
    return HttpResponse("Hello World! - hello page")
    #return render(request, 'home.html') #, {'documents': documents, 'form': form})

def home(request):
    
    return render(request, 'home.html') #, {'documents': documents, 'form': form})
    #return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def python_code(request):
    
    return render(request, 'python_code.html') #, {'documents': documents, 'form': form})
    #return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def frequent_questions(request):
    
    return render(request, 'frequent_questions.html') #, {'documents': documents, 'form': form})
    #return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def email(request):
    
    form_class = ContactForm

    # new logic!
    if request.method == 'POST':
        form = form_class(data=request.POST)

        if form.is_valid():
            contact_name = request.POST.get(
                'contact_name'
            , '')
            contact_email = request.POST.get(
                'contact_email'
            , '')
            form_content = request.POST.get('content', '')

            # Email the profile with the 
            # contact information
            template = loader.get_template('contact_template.txt')
            #context = Context({
            #    'contact_name': contact_name,
            #    'contact_email': contact_email,
            #    'form_content': form_content,
            #}) 
            #content = template.render(context)

            email = EmailMessage(
                "New Email",
                "From: "+contact_name+"\n\n"+
                "Email: "+contact_email+"\n\n"+
                form_content,
                'cat@alumni.ubc.ca',
                ['catubc@catubc.webfactional.com'],
                reply_to=['catubc@catubc.webfactional.com']
            )
            #email = EmailMessage('Hello', 'Body goes here', 'from@example.com',
            #['to1@example.com', 'to2@example.com'], ['bcc@example.com'],
            #reply_to=['another@example.com'], headers={'Message-ID': 'foo'})
            
            email.send()

            #return redirect('http://www.spikesortingtest.com/email')

    return render(request, 'email.html', {
        'form': form_class,
    })

def send_email(request):
    msg = EmailMessage('Request Callback',
                       'Here is the message.', to=['charl@byteorbit.com'])
    msg.send()
    return HttpResponseRedirect('/')
    
    #return render(request, 'email.html') #, {'documents': documents, 'form': form})

def docs(request):    
    
    return render(request, 'docs.html') #, {'documents': documents, 'form': form})
    #return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def metrics(request):
    
    #Call routine that will compute .jpgs for metrics webpage
    Compute_global_metrics()
        
    return render(request, 'metrics.html') #, {'documents': documents, 'form': form})
    #return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def delete_matrix():
    documents = Document.objects.all()
    for document in documents:
        document.delete()


def list_dir(request):
    # pdb.set_trace()

    # my_csv_list = MyCSVImporterModel(source="/home/cornelis/PycharmProjects/minimal-django-file-upload-example/src/for_django_1-8/myproject/media/documents/2016/01/16/silico_0_truth_filtered_ptps.csv")
    # row, first_line = my_csv_list.cleaned_data[0]
    # first_line['age']

    # ---------- from minimal-django-file-upload-example
    newdoc = None
    # Handle file upload
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
           # pdb.set_trace()
            #todo: Find out how to delete previousy uploaded files
            #delete_matrix()
            newdoc = Document(docfile=request.FILES['docfile'])
            newdoc.save()

            # Redirect to the document list after POST
            # TODO: Find out what this does exactly (comment out still "works")
            # return HttpResponseRedirect(reverse('myproject.myapp.views.list'))
    else:
        form = DocumentForm()  # A empty, unbound form

    documents = Document.objects.all()
    csvfile = None
    if len(documents) > 0:
        if newdoc is None:
            x=1
            #plot = figure()
            #script, div = components(plot, CDN)
        else:
            csvfile = newdoc.docfile.path

            return plotResult(request, csvfile)

    return render(request, 'list_dir.html', {'documents': documents, 'form': form})

def datasets(request):

    # ---------- from minimal-django-file-upload-example
    newdoc = None
    # Handle file upload
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            newdoc = Document(docfile=request.FILES['docfile'])
            newdoc.save()
    else:
        form = DocumentForm()  # A empty, unbound form

    documents = Document.objects.all()
    csvfile = None
    if len(documents) > 0:
        if newdoc is None:
            x=1
            #plot = figure()
            #script, div = components(plot, CDN)
        else:
            csvfile = newdoc.docfile.path

            return plotResult(request, csvfile)

    return render(request, 'datasets.html', {'documents': documents, 'form': form})

def upload_current(request):

    # ---------- from minimal-django-file-upload-example
    newdoc = None
    # Handle file upload
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            newdoc = Document(docfile=request.FILES['docfile'])
            newdoc.save()
    else:
        form = DocumentForm()  # A empty, unbound form

    documents = Document.objects.all()
    csvfile = None
    if len(documents) > 0:
        if newdoc is None:
            pass
            #x=1
            #plot = figure()
            #script, div = components(plot, CDN)
        else:
            csvfile = newdoc.docfile.path
            
            #Invert the dictionary request.POST to reverse lookup
            inv_map = {v: k for k, v in request.POST.items()}
            
            truth_file = inv_map['Upload']
            #return HttpResponse(truth_file)
            
            #return plotResult(request, csvfile, truth_file)

            random_file = plotResult(request, csvfile, truth_file) #Save .jpg to directory
            
            return render(request, 'output.html', {"variable_int" : random_file})
            
            #from django.shortcuts import redirect
            #return redirect('http://www.spikesortingtest.com/output')

    return render(request, 'upload.html', {'uploads': documents, 'form': form})


def upload(request):

    # ---------- from minimal-django-file-upload-example
    newdoc = None
    # Handle file upload
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            newdoc = Document(docfile=request.FILES['docfile'])
            newdoc.save()
    else:
        form = DocumentForm()  # A empty, unbound form

    documents = Document.objects.all()
    csvfile = None
    if len(documents) > 0:
        if newdoc is None:
            pass

        else:
            csvfile = newdoc.docfile.path
            
            #Invert the dictionary request.POST to reverse lookup
            inv_map = {v: k for k, v in request.POST.items()}
            
            truth_file = inv_map['Upload']
            #return HttpResponse(truth_file)
            
            #return plotResult(request, csvfile, truth_file)

            #plotResult(csvfile, truth_file) #Check the sort and generate outputs
            
            #import threading
            #t = threading.Thread(target=plotResult(csvfile, truth_file))
            #t.setDaemon(True)
            #t.start()
            ##return HttpResponse()
            
            random_file = csvfile.replace('/home/catubc/webapps/my_django_app/myproject/media/uploads/','')
            print "...random_file: ", random_file
            makeEmptyResults(random_file, truth_file)

            #random_file = plotResult(request, csvfile, truth_file) #Save .jpg to directory
            
            return render(request, 'output.html', {"uploaded_file" : random_file, "selected_file" : truth_file})
            
            #from django.shortcuts import redirect
            #return redirect('http://www.spikesortingtest.com/output')

    return render(request, 'upload.html', {'uploads': documents, 'form': form})


def makeEmptyResults(random_file, truth_file):
    
    from shutil import copyfile
    black_jpg_filename = '/home/catubc/webapps/static_media/images/unprocessed.jpg'
    saved_output_jpg = '/home/catubc/webapps/static_media/saved_sort_images/' + random_file+'.jpg'
    copyfile(black_jpg_filename,saved_output_jpg)
    
    #Save purity and completeness .txt files
    data_empty = np.zeros(2, dtype=np.int16)
    np.savetxt('/home/catubc/webapps/static_media/saved_sort_images/purity_' + random_file+'.txt', data_empty)
    np.savetxt('/home/catubc/webapps/static_media/saved_sort_images/completeness_' + random_file+'.txt', data_empty)


def process_sort(request):
    
    print "...type: ", type(request.POST.get('your_filename'))
    #if type(request.POST.get('your_filename')) is None:
    if (request.POST.get('your_filename') is None) or (len(request.POST.get('your_filename')) == 0):
        return render(request, 'output.html', {"uploaded_file" : '', "selected_file" : ''})

    else:
        your_filename = '/home/catubc/webapps/my_django_app/myproject/media/uploads/'+request.POST.get('your_filename')
        selected_simulation = request.POST.get('selected_simulation')
        print "....input files: ", your_filename, selected_simulation

        #Check to see if data already processed: 
        purity_data = np.loadtxt('/home/catubc/webapps/static_media/saved_sort_images/purity_' + request.POST.get('your_filename')+'.txt')
        print purity_data
        if purity_data[0]==0.0: 
            print "... computing sorting metrics..."
            plotResult(your_filename, selected_simulation) #Check the sort and generate outputs

        return render(request, 'output.html', {"uploaded_file" : request.POST.get('your_filename'), "selected_file" : selected_simulation})
    

#def upload_backup(request):

    ## ---------- from minimal-django-file-upload-example
    #newdoc = None
    ## Handle file upload
    #if request.method == 'POST':
        #form = DocumentForm(request.POST, request.FILES)
        #if form.is_valid():
            #newdoc = Document(docfile=request.FILES['docfile'])
            #newdoc.save()
    #else:
        #form = DocumentForm()  # A empty, unbound form

    #documents = Document.objects.all()
    #csvfile = None
    #if len(documents) > 0:
        #if newdoc is None:
            #pass
            ##x=1
            ##plot = figure()
            ##script, div = components(plot, CDN)
        #else:
            #csvfile = newdoc.docfile.path
            
            ##Invert the dictionary request.POST to reverse lookup
            #inv_map = {v: k for k, v in request.POST.items()}
            
            #truth_file = inv_map['Upload']
            ##return HttpResponse(truth_file)

            #return plotResult(request, csvfile, truth_file)
            
            ##from django.shortcuts import redirect
            ##return redirect('http://www.spikesortingtest.com/output')

    #return render(request, 'upload.html', {'uploads': documents, 'form': form})

def output(request):

    return render(request, 'output.html')


def plotResult(sorted_dir, truth_file):

    sorted_file = str(truth_file)
    data_dir = settings.STATICFILES_DIRS[0]+'/'+sorted_file+'/'

    #********************* READ TSF FILES *******************
    tsf_name = data_dir + sorted_file  + '_truth_filtered'  #NB: Not loading channel data, only header info

    f = open(tsf_name + '.tsf', "rb")
    print "Loading ", tsf_name + '.tsf'
    tsf = Tsf_file(f, '')  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
    tsf.tsf_name = tsf_name
    tsf.sorted_file = sorted_file+ '_truth_filtered'
    f.close()

    #else: 
        #tsf_name = data_dir + '/'+ sorted_file + '_truth_filtered'
        #class Tsf: 
            #pass
        #tsf = Tsf() #MAKE DUMMY OBJECT

    ##*************** READ GROUND TRUTH FILE : SORT1 ************************
    ptcs_flag = 1
    if (sorted_file == "ViSAPy_0"):
        Sort1 = dirk_2.Loadcsv_vispy(sorted_file, data_dir, tsf)
        #file_name = '/media/cat/12TB/spikesortingtest.com/espen/ViSAPy_ground_truth.gdf'

    #Allen Inst biophys simulations
    else: 
        Sort1 = dirk_2.Loadptcs(sorted_file+'_truth_filtered', data_dir, ptcs_flag, save_timestamps=False)
        Sort1.name=sorted_file+'_truth_filtered'
        Sort1.filename=sorted_file+'_truth_filtered'
        Sort1.rawfile= tsf_name
        Sort1.directory=data_dir
        Sort1.data_dir = data_dir
    #Vispy biophys simulations; have own .gdf spike time format.
        
    ##************** READ 2ND SORT DATA **********
    fname_save = sorted_dir
    Sort2 = dirk_2.Loadcsv(sorted_dir, Sort1, sorted_file + '/')
    Sort2.directory= sorted_dir
    Sort2.name = sorted_file
    Sort2.filename=sorted_file
    Sort2.chanpos=[999]

    Sort2.n_spikes=0
    for i in range(len(Sort2.units)):
        Sort2.n_spikes+= len(Sort2.units[i])

    #Run the compare_sorts algorithm on the data if not already done
    dirk_2.Compare_sorts(Sort1, Sort2) #Run the sorter;save purity and completeness inside Sort2.

    # todo: somehow store image in a string buffer and send buffer in http response to browser

    return Plot_Composite_metrics(data_dir, sorted_file, fname_save, Sort2)

def simple_chart(request):
    from bokeh.plotting import figure
    from bokeh.resources import CDN
    from bokeh.embed import file_html

    plot = figure()
    plot.circle([1,2], [3,4])

    html = file_html(plot, CDN, "my plot")
    return HttpResponse(html)

