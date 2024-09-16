# models.py
from peewee import *

databaseFile='/data/diriano/sugarcaneMapperIDs/sugarcaneSequences.db'
# Creating the database
db = SqliteDatabase(databaseFile, pragmas={
    'foreign_keys': 1,  # Enforce foreign-key constraints
    'ignore_check_constraints': 0,  # Enforce CHECK constraints
    'journal_mode':'WAL'  # Use WAL mode for better concurrency,
})

class BaseModel(Model):
    class Meta:
        database = db

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

class SequenceSet(BaseModel):
    seID = PrimaryKeyField()
    nameSet = CharField(unique=True)

    class Meta:
        db_table = 'sequenceSet'

class Sequence2Set(BaseModel):
    sequenceID = ForeignKeyField(Sequence, column_name='sequenceID', backref='sequence2sets', on_delete='CASCADE')
    seID = ForeignKeyField(SequenceSet, column_name='seID', backref='sequence2sets', on_delete='CASCADE')

    class Meta:
        db_table = 'sequence2Set'
        primary_key = CompositeKey('sequenceID', 'seID')

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