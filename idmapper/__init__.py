import os

from flask import Flask, render_template, g
from markupsafe import escape
from idmapper.models import db, Sequence, panTranscriptomeGroup

def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # Connect to the database at the start of each request
    @app.before_request
    def before_request():
        g.db = db
        if g.db.is_closed():
            g.db.connect()

    # Close the database after the request
    @app.teardown_request
    def teardown_request(exception=None):
        if not g.db.is_closed():
            g.db.close()
            
    # a simple page that says hello
    @app.route('/')
    def index():
        return 'Welcome to the IDMapper for Conekt Grasses'

    @app.route('/mapid/<id>')
    def mapid(id=None):
        try:
            sequence = Sequence.get(Sequence.sequenceIdentifier == id)
            return render_template('idmap.html', identifier=sequence)
        except Sequence.DoesNotExist:
            return f"No sequence found with ID: {id}", 404

    return app
