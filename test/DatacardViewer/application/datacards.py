import os
from datetime import timedelta
from application import app
from flask import make_response, send_file, request, Response, url_for, jsonify, send_from_directory, current_app
from json import JSONEncoder
from werkzeug import secure_filename
from functools import update_wrapper
import paramiko

# This is the path to the upload directory
app.config['UPLOAD_FOLDER'] = 'application/uploads/'
# These are the extension that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['txt','root'])
# Maximum size
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

# def get_file_size(file):
#     file.seek(0, 2)  # Seek to the end of the file
#     size = file.tell()  # Get the position of EOF
#     file.seek(0)  # Reset the file position to the beginning
#     return size

def _handleUpload(files):
    if not files:
       return None
    filenames = []
    saved_files_urls = []
    for key, file in files.iteritems():
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            saved_files_urls.append(url_for('uploaded_file', filename=filename))
            filenames.append("%s" % (file.filename))
    return filenames

@app.route('/datacards', methods=['POST'])
def upload_datacards():
    try:
        files = request.files
        uploaded_files = _handleUpload(files)
        return jsonify({'files': uploaded_files})
    except:
        raise
        return jsonify({'status': 'error'})

'''
#TODO Retrieve list of uploaded datacards [TO SHOW AFTER UPLOADS]
@app.route('/datacards', methods = ['GET'])
def get_datacards_list():
    filenames = []
    for filename in os.listdir(app.config['UPLOAD_FOLDER']):
        filenames.append(filename)
    return jsonify({'files': filenames})
'''
#Retrieve list of datacards 
@app.route('/datacards/list', methods = ['GET'])
def get_datacards_list():
    filenames = []
    for filename in os.listdir(app.config['UPLOAD_FOLDER']):
        filename, ext = os.path.splitext(filename)
        if ext == '.txt':
            filenames.append(filename+ext)
    return jsonify({'files': filenames})

from DatacardParser import *
@app.route('/datacards/<filename>')
def uploaded_file(filename):
    from optparse import OptionParser

    parser = OptionParser()
    addDatacardParserOptions(parser)
    (options, args) = parser.parse_args()
    options.out = "tmp.root"
    options.fileName = filename
    options.allowNoSignal = True
    fileToRead = open(app.config['UPLOAD_FOLDER']+filename, 'r')
    datacard = parseCard(fileToRead, options)

    jsonString = JSONEncoder().encode({
        "filename": filename,
        "observation": datacard.obs,
        "binsProcessesRates": datacard.exp,
        "nuisances": datacard.systs,
        "shapeMap": datacard.shapeMap
    })
    return jsonString


#-- partial content streaming
# from flask import request, send_file, Response

@app.after_request
def after_request(response):
    response.headers.add('Accept-Ranges', 'bytes')
    return response

def send_file_partial(path, **kwargs):
    """
    Simple wrapper around send_file which handles HTTP 206 Partial Content
    (byte ranges)
    TODO: handle all send_file args, mirror send_file's error handling
    (if it has any)
    """
    range_header = request.headers.get('Range', None)
    if not range_header:
        return make_response(open(path).read())

    size = os.path.getsize(path)
    byte1, byte2 = 0, None
    
    m = re.search('(\d+)-(\d*)', range_header)
    g = m.groups()
    
    if g[0]:
        byte1 = int(g[0])

    if g[1]:
        byte2 = int(g[1])

    length = size - byte1
    if byte2 is not None:
        length = byte2 - byte1
    
    data = None
    with open(path, 'rb') as f:
        f.seek(byte1)
        data = f.read(length)

    rv = Response(data,
        206,
        mimetype='application/root+root.exe',
        direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(byte1, byte1 + length - 1, size))

    return rv

@app.route('/datacards/files/<filename>', methods = ['GET','OPTIONS'])
def get_file(filename):
    return send_file_partial('application/uploads/'+filename, cache_timeout=0 )