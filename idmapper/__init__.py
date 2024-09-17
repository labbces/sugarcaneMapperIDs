import os

from flask import Flask, render_template, g
from markupsafe import escape
from idmapper.models import db, Sequence, panTranscriptomeGroup, Sequence2Set

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
            # Get the sequence by its identifier
            # Selecting specific fields from the Sequence model
            proteinSequence = Sequence.select(Sequence.sequenceIdentifier).where(
                (Sequence.sequenceIdentifier == id) & 
                (Sequence.sequenceVersion == 1) & 
                (Sequence.sequenceClass == 'protein')
            ).get()
        
            # Check if the sequence belongs to a panTranscriptomeGroup
            pan_group = panTranscriptomeGroup.get_or_none(panTranscriptomeGroup.sequenceID == proteinSequence.ID)
            
            # Check if the sequence belongs to a Sequence2Set
            sequence_set = Sequence2Set.get_or_none(Sequence2Set.sequenceID == proteinSequence.ID)
            
            return render_template(
                'idmap.html', 
                identifier=id, 
                sequence=proteinSequence, 
                pan_group=pan_group, 
                sequence_set=sequence_set
            )
        except Sequence.DoesNotExist:
            return render_template('idmap.html', identifier=id, sequence=None)

    return app
