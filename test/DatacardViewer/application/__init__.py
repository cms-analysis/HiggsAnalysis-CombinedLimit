from flask import Flask

# Initialize the Flask application
app = Flask(__name__)

from application import views
from application import datacards