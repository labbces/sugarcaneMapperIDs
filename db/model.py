# models.py
from peewee import *

# Creating the database
db = SqliteDatabase('sugarcaneSequences.db', pragmas={
    'foreign_keys': 1,  # Enforce foreign-key constraints
    'ignore_check_constraints': 0,  # Enforce CHECK constraints
})

class BaseModel(Model):
    class Meta:
        database = db

# The Sequence table stores DNA (gene, transcript, and CDS) and protein sequences
class Sequence(BaseModel):
    ID = PrimaryKeyField()
    sequenceIdentifier = CharField(index=True)
    sequenceVersion = IntegerField(default=1)
    sequenceType = CharField(max_length=255, default='DNA', constraints=[SQL('CHECK( sequenceType IN ("DNA", "Amino acid"))')])
    sequenceClass = CharField(max_length=255, null=False, default='gene', constraints=[SQL('CHECK( sequenceClass IN ("gene", "transcript", "CDS", "protein"))')])
    sequenceLength = IntegerField(null=False)
    sequence = TextField(null=True)

    class Meta:
        db_table = 'seq'
        constraints = [SQL('CONSTRAINT uniqueSequenceIdVersionType UNIQUE(sequenceIdentifier,sequenceVersion,sequenceClass) ON CONFLICT ABORT')]

# The SequenceSet table stores the name of the sets of sequences
class SequenceSet(BaseModel):
    seID = PrimaryKeyField()
    nameSet = CharField(unique=True)

    class Meta:
        db_table = 'sequenceSet'

# Link the sequence to the SequenceSets
class Sequence2Set(BaseModel):
    sequenceID = ForeignKeyField(Sequence, column_name='sequenceID', backref='sequence2sets', on_delete='CASCADE')
    seID = ForeignKeyField(SequenceSet, column_name='seID', backref='sequence2sets', on_delete='CASCADE')

    class Meta:
        db_table = 'sequence2Set'
        primary_key = CompositeKey('sequenceID', 'seID')

# The panTranscriptomeGroup table stores the groups of sequences
class panTranscriptomeGroup(BaseModel):
    panTranscriptomeGroupID = PrimaryKeyField()
    sequenceID = ForeignKeyField(Sequence, column_name='sequenceID', backref='pantranscriptome_groups', on_delete='CASCADE')
    groupID = CharField()
    representative = BooleanField(default=False)

    class Meta:
        db_table = 'panTranscriptomeGroup'
        constraints = [SQL('CONSTRAINT uniqueSequenceIdGroupId UNIQUE(sequenceID,groupID)')]

# List of models
tables = [Sequence, SequenceSet, Sequence2Set, panTranscriptomeGroup]

db.connect()

try:
    for table in tables:
        if not table.table_exists():
            db.create_tables([table])
            print(f"Table '{table.__name__}' created succesfully!")
        else:
            print(f"Table '{table.__name__}' already exists!")
except OperationalError as e:
    print(f"Error occurred while creating the tables: {e}")
except Exception as e:
    print(f"Error occurred: {e}")
finally:
    db.close()
