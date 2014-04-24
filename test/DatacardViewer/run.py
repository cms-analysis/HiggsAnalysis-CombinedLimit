#!flask/bin/python
from application import app

app.run(
		debug = True, 
		host = '137.138.216.35', 
		port = 5000
	)