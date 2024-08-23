PRAGMA foreign_keys = ON;
PRAGMA ignore_check_constraints = OFF;
DROP TABLE IF EXISTS sequence;
DROP TABLE IF EXISTS sequenceSet;
DROP TABLE IF EXISTS sequence2sequenceSet;

CREATE TABLE seq(
  ID INTEGER PRIMARY KEY AUTOINCREMENT,
  sequenceIdentifier VARCHAR(255) NOT NULL,
  sequenceVersion INTEGER NOT NULL,
  sequenceType VARCHAR(255) CHECK ( sequenceType IN ('DNA', 'Protein') ) NOT NULL DEFAULT 'DNA',
  sequence TEXT NOT NULL,
  CONSTRAINT uniqueSequenceIdVersionType UNIQUE(sequenceIdentifier,sequenceVersion,sequenceType) ON CONFLICT ABORT
);

CREATE TABLE sequenceSet (
  seID INTEGER PRIMARY KEY AUTOINCREMENT,
  nameSet VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE sequence2set (
  sequenceID INTEGER NOT NULL,
  seID INTEGER NOT NULL,
  PRIMARY KEY(sequenceID,seID),
  FOREIGN KEY(sequenceID) REFERENCES seq(ID),
  FOREIGN KEY (seID) REFERENCES sequenceSet(seID)
);

CREATE TABLE panTranscriptomeGroup (
  ID INTEGER PRIMARY KEY AUTOINCREMENT,
  sequenceID INTEGER NOT NULL,
  groupID VARCHAR(255) NOT NULL,
  representative BOOLEAN NOT NULL DEFAULT 0,
  FOREIGN KEY(sequenceID) REFERENCES seq(ID),
  UNIQUE (sequenceID,groupID)
);

INSERT INTO seq (sequenceIdentifier, sequenceVersion, sequenceType, sequence) VALUES ('ppppp2', 1, 'DNA', 'acacacacacacacac');
INSERT INTO seq (sequenceIdentifier, sequenceVersion, sequenceType, sequence) VALUES ('ppppp2', 2, 'DNA', 'atacacacacacacac');
INSERT INTO seq (sequenceIdentifier, sequenceVersion, sequenceType, sequence) VALUES ('ppppp1', 1, 'DNA', 'atgcacacacacacac');

INSERT INTO sequenceSet (nameSet) VALUES ('set1');
INSERT INTO sequenceSet (nameSet) VALUES ('set2');

INSERT INTO sequence2set (sequenceID, seID) VALUES (1, 1);
INSERT INTO sequence2set (sequenceID, seID) VALUES (2, 2);
INSERT INTO sequence2set (sequenceID, seID) VALUES (3, 1);

INSERT INTO panTranscriptomeGroup (sequenceID, groupID, representative) VALUES (1, 'group1', 1);
INSERT INTO panTranscriptomeGroup (sequenceID, groupID, representative) VALUES (3, 'group1', 0);
INSERT INTO panTranscriptomeGroup (sequenceID, groupID, representative) VALUES (2, 'group2', 1);