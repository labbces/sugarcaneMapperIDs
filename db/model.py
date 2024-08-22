from  peewee import *

#creating the database
# this db will store all sugarcane sequences, their ids and their sequence, as well the the orthogroups confomration and orthogroup representative sequences that we use in conekt
db = SqliteDatabase('sugarcaneSequences.db', pragmas={
    'foreign_keys': 1,  # Enforce foreign-key constraints
    'ignore_check_constraints': 0,  # Enforce CHECK constraints
    })

class BaseModel(Model):
    class Meta:
        database = db

#The sequence table store DNA (gene, transcript and CDS) and protein sequences
class Sequence(BaseModel):
# Create a class based on the following table definition:
#   CREATE TABLE sequence(
#   ID INTEGER PRIMARY KEY AUTOINCREMENT,
#   sequenceIdentifier VARCHAR(255) NOT NULL,
#   sequenceVersion INTEGER NOT NULL,
#   sequenceType VARCHAR(255) CHECK ( sequenceType IN ('DNA', 'Protein') ) NOT NULL DEFAULT 'DNA',
#   sequence TEXT NOT NULL,
#   CONSTRAINT uniqueSequenceIdVersionType,
#   UNIQUE(sequenceIdentifier,sequenceVersion,sequenceType) ON CONFLICT ABORT
#   );

#   CREATE INDEX idx_sequenceIdentifier ON sequence(sequenceIdentifier);

    ID = PrimaryKeyField()
    sequenceIdentifier = CharField(index=True)
    sequenceVersion = IntegerField(default=1)
    sequenceType = CharField(max_length=255, default='DNA', constraints=[SQL('CHECK( sequenceType IN ("DNA", "Protein"))')])
    sequence = TextField()

    class Meta:
        db_table = 'sequence'
        constraints = [SQL('CONSTRAINT uniqueSequenceIdVersionType UNIQUE(sequenceIdentifier,sequenceVersion,sequenceType) ON CONFLICT ABORT')]

#The SequenceSet table store the name of the sets of sequences that we use in conekt, genotypes, pantranscriptome, eventually genome versions
class SequenceSet(BaseModel):
# Create a class based on the following table definition:
#   CREATE TABLE sequenceSet (
#   ID INTEGER PRIMARY KEY AUTOINCREMENT,
#   nameSet VARCHAR(255) NOT NULL UNIQUE,
#   );
# CREATE INDEX idx_nameSet ON sequenceSet(nameSet);
    ID = PrimaryKeyField()
    nameSet = CharField(unique=True)

    class Meta:
        db_table = 'sequenceSet'

#Link the sequence to the SequenceSets
class Sequence2sequenceSet(BaseModel):
# Create a class based on the following table definition:
# CREATE TABLE sequence2sequenceSet (
#   sequenceID INTEGER NOT NULL,
#   sequenceSetID INTEGER NOT NULL
#   PRIMARY KEY(sequenceID,sequenceSetID),
#   FOREIGN KEY(sequenceID) REFERENCES sequence(ID),
#   FOREIGN KEY (sequenceSetID) REFERENCES sequenceSet(ID)
# );
    sequenceID = ForeignKeyField(Sequence, field='ID')
    sequenceSetID = ForeignKeyField(SequenceSet, field='ID')

    class Meta:
        db_table = 'sequence2sequenceSet'
        primary_key = CompositeKey('sequenceID', 'sequenceSetID')

#List of models
tables=[Sequence, SequenceSet, Sequence2sequenceSet]

db.connect()

try:
    for table in tables:
        if not table.table_exists():
            db.create_tables([table])
            print(f"Table '{table.__name__}' created succesfully!")
        else:
            print(f"Table '{table.__name__}' already exists!")
except OperationalError as e:
    print(f"Error occurred while creating the table: {e}")
except Exception as e:
    print(f"Error occurred: {e}")
finally:
    db.close()


