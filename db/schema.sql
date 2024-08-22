DROP TABLE IF EXISTS sequence;
DROP TABLE IF EXISTS sequenceSet;
DROP TABLE IF EXISTS sequence2sequenceSet;

CREATE TABLE sequence(
  ID INTEGER PRIMARY KEY AUTOINCREMENT,
  sequenceIdentifier VARCHAR(255) NOT NULL,
  sequenceVersion INTEGER NOT NULL,
  sequenceType VARCHAR(255) CHECK ( sequenceType IN ('DNA', 'Protein') ) NOT NULL DEFAULT 'DNA',
  sequence TEXT NOT NULL,
  CONSTRAINT uniqueSequenceIdVersionType,
  UNIQUE(sequenceIdentifier,sequenceVersion,sequenceType) ON CONFLICT ABORT
);

CREATE INDEX idx_sequenceIdentifier ON sequence(sequenceIdentifier);

CREATE TABLE sequenceSet (
  ID INTEGER PRIMARY KEY AUTOINCREMENT,
  nameSet VARCHAR(255) NOT NULL UNIQUE,
);

CREATE INDEX idx_nameSet ON sequenceSet(nameSet);

CREATE TABLE sequence2sequenceSet (
  sequenceID INTEGER NOT NULL,
  sequenceSetID INTEGER NOT NULL
  PRIMARY KEY(sequenceID,sequenceSetID),
  FOREIGN KEY(sequenceID) REFERENCES sequence(ID),
  FOREIGN KEY (sequenceSetID) REFERENCES sequenceSet(ID)
);

CREATE TABLE panTranscriptome (
  ID INTEGER PRIMARY KEY AUTOINCREMENT,
  sequenceID INTEGER NOT NULL,
  groupID VARCHAR(255) NOT NULL,
  FOREIGN KEY(sequenceSetID) REFERENCES sequence(ID),
  PRIMARY KEY(sequenceID,groupID)
);