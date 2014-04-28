#!flask/bin/python
from application import app
import urllib2

app.run(
		debug = True, 
		host = urllib2.urlopen("http://myip.dnsdynamic.org/").read(), 
		port = 5000
	)
