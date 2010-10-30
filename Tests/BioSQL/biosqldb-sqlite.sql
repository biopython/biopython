--  BioSQL database schema for SQLite.
-- 
--  This file is part of BioSQL.
--
--  BioSQL is free software: you can redistribute it and/or modify it
--  under the terms of the GNU Lesser General Public License as
--  published by the Free Software Foundation, either version 3 of the
--  License, or (at your option) any later version.
--
--  BioSQL is distributed in the hope that it will be useful,
--  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  GNU Lesser General Public License for more details.
--
--  You should have received a copy of the GNU Lesser General Public License
--  along with BioSQL. If not, see <http://www.gnu.org/licenses/>.
--
-- ========================================================================
--
-- See MySQL database schema and BioSQL website for table documentation.
-- This contains notes specific to SQLite
-- 
-- A note about Primary Keys in SQLite
-- SQLite automatically creates a ROWID for each row of a table.
--   Using this ROWID as the primary key is faster than using a 
--   user-defined primary key. By declaring a column as an 
--   INTEGER PRIMARY KEY, you are actually creating an alias to the
--   ROWID and get the associated speed benefits.  The ROWID effectively
--   "autoincrements"; however, it can reuse ROWIDs of deleted rows. 
--   To avoid reusing old ROWIDs would require adding the AUTOINCREMENT
--   keyword, which also reduces the performance.
--   ( see http://www.sqlite.org/autoinc.html) 

CREATE TABLE biodatabase (
  	biodatabase_id 	INTEGER PRIMARY KEY,
  	name           	VARCHAR(128) NOT NULL,
	authority	VARCHAR(128),
	description	TEXT,
  	UNIQUE (name)
);

CREATE INDEX db_auth on biodatabase(authority);

CREATE TABLE taxon (
       taxon_id		INTEGER PRIMARY KEY,
       ncbi_taxon_id 	INT(10),
       parent_taxon_id	INT(10) ,
       node_rank	VARCHAR(32),
       genetic_code	TINYINT ,
       mito_genetic_code TINYINT ,
       left_value	INT(10) ,
       right_value	INT(10) ,
       UNIQUE (ncbi_taxon_id),
       UNIQUE (left_value),
       UNIQUE (right_value)
);

CREATE INDEX taxparent ON taxon(parent_taxon_id);

CREATE TABLE taxon_name (
       taxon_id		INTEGER,
       name		VARCHAR(255)  NOT NULL,
       name_class	VARCHAR(32)  NOT NULL,
       UNIQUE (taxon_id,name,name_class)
);

CREATE INDEX taxnametaxonid ON taxon_name(taxon_id);
CREATE INDEX taxnamename    ON taxon_name(name);

CREATE TABLE ontology (
       	ontology_id        INTEGER PRIMARY KEY,
       	name	   	   VARCHAR(32)  NOT NULL,
       	definition	   TEXT,
	UNIQUE (name)
);

CREATE TABLE term (
       	term_id   INTEGER PRIMARY KEY,
       	name	   	   VARCHAR(255)  NOT NULL,
       	definition	   TEXT,
	identifier	   VARCHAR(40) ,
	is_obsolete	   CHAR(1),
	ontology_id	   INTEGER,
	UNIQUE (identifier),
        UNIQUE (name,ontology_id,is_obsolete)
);

CREATE INDEX term_ont ON term(ontology_id);

CREATE TABLE term_synonym (
       synonym		  VARCHAR(255)  NOT NULL,
       term_id		  INTEGER,
       PRIMARY KEY (term_id,synonym)
);

CREATE TABLE term_dbxref (
       	term_id	          INTEGER,
       	dbxref_id         INTEGER,
	rank		  SMALLINT,
	PRIMARY KEY (term_id, dbxref_id)
);

CREATE INDEX trmdbxref_dbxrefid ON term_dbxref(dbxref_id);

CREATE TABLE term_relationship (
        term_relationship_id INTEGER PRIMARY KEY,
       	subject_term_id	INTEGER,
       	predicate_term_id    INTEGER,
       	object_term_id       INTEGER,
	ontology_id	INTEGER,
	UNIQUE (subject_term_id,predicate_term_id,object_term_id,ontology_id)
);

CREATE INDEX trmrel_predicateid ON term_relationship(predicate_term_id);
CREATE INDEX trmrel_objectid ON term_relationship(object_term_id);
CREATE INDEX trmrel_ontid ON term_relationship(ontology_id);

CREATE TABLE term_relationship_term (
        term_relationship_id INTEGER PRIMARY KEY,
        term_id              INTEGER,
        UNIQUE ( term_id ) 
);

CREATE TABLE term_path (
        term_path_id         INTEGER PRIMARY KEY,
       	subject_term_id	     INTEGER,
       	predicate_term_id    INTEGER,
       	object_term_id       INTEGER,
	ontology_id          INTEGER,
	distance	     INT(10) ,
	UNIQUE (subject_term_id,predicate_term_id,object_term_id,ontology_id,distance)
);

CREATE INDEX trmpath_predicateid ON term_path(predicate_term_id);
CREATE INDEX trmpath_objectid ON term_path(object_term_id);
CREATE INDEX trmpath_ontid ON term_path(ontology_id);

CREATE TABLE bioentry (
	bioentry_id	    INTEGER PRIMARY KEY,
  	biodatabase_id  INTEGER,
  	taxon_id     	INT(10) ,
  	name		VARCHAR(40) NOT NULL,
  	accession    	VARCHAR(128)  NOT NULL,
  	identifier   	VARCHAR(40) ,
	division	VARCHAR(6),
  	description  	TEXT,
  	version 	SMALLINT  NOT NULL, 
  	UNIQUE (accession,biodatabase_id,version),
 	UNIQUE (identifier, biodatabase_id)
);

CREATE INDEX bioentry_name ON bioentry(name);
CREATE INDEX bioentry_db   ON bioentry(biodatabase_id);
CREATE INDEX bioentry_tax  ON bioentry(taxon_id);

CREATE TABLE bioentry_relationship (
        bioentry_relationship_id INTEGER PRIMARY KEY,
        object_bioentry_id 	 INTEGER,
   	subject_bioentry_id 	 INTEGER,
   	term_id 		 INTEGER,
   	rank 			 INT(5),
	UNIQUE (object_bioentry_id,subject_bioentry_id,term_id)
);

CREATE INDEX bioentryrel_trm   ON bioentry_relationship(term_id);
CREATE INDEX bioentryrel_child ON bioentry_relationship(subject_bioentry_id);

CREATE TABLE bioentry_path (
   	object_bioentry_id 	INTEGER PRIMARY KEY,
   	subject_bioentry_id 	INTEGER,
   	term_id 		INTEGER,
	distance	     	INT(10) ,
	UNIQUE (object_bioentry_id,subject_bioentry_id,term_id,distance)
);

CREATE INDEX bioentrypath_trm   ON bioentry_path(term_id);
CREATE INDEX bioentrypath_child ON bioentry_path(subject_bioentry_id);

CREATE TABLE biosequence (
  	bioentry_id     INTEGER PRIMARY KEY,
  	version     	SMALLINT, 
  	length      	INT(10),
  	alphabet        VARCHAR(10),
  	seq 		LONGTEXT
);

CREATE TABLE dbxref (
        dbxref_id	INTEGER PRIMARY KEY,
        dbname          VARCHAR(40)  NOT NULL,
        accession       VARCHAR(128)  NOT NULL,
	version		SMALLINT  NOT NULL,
        UNIQUE(accession, dbname, version)
);

CREATE INDEX dbxref_db  ON dbxref(dbname);

CREATE TABLE dbxref_qualifier_value (
       	dbxref_id 		INTEGER,
       	term_id 		INTEGER,
  	rank  		   	SMALLINT NOT NULL DEFAULT 0,
       	value			TEXT,
	PRIMARY KEY (dbxref_id,term_id,rank)
);

CREATE INDEX dbxrefqual_dbx ON dbxref_qualifier_value(dbxref_id);
CREATE INDEX dbxrefqual_trm ON dbxref_qualifier_value(term_id);

CREATE TABLE bioentry_dbxref ( 
       	bioentry_id        INTEGER,
       	dbxref_id          INTEGER,
  	rank  		   SMALLINT,
	PRIMARY KEY (bioentry_id,dbxref_id)
);

CREATE INDEX dblink_dbx  ON bioentry_dbxref(dbxref_id);

CREATE TABLE reference (
  	reference_id       INTEGER PRIMARY KEY,
	dbxref_id	   INT(10) ,
  	location 	   TEXT NOT NULL,
  	title    	   TEXT,
  	authors  	   TEXT,
  	crc	   	   VARCHAR(32),
	UNIQUE (dbxref_id),
	UNIQUE (crc)
);

CREATE TABLE bioentry_reference (
  	bioentry_id 	INTEGER,
  	reference_id 	INTEGER,
  	start_pos	INT(10),
  	end_pos	  	INT(10),
  	rank  		SMALLINT NOT NULL DEFAULT 0,
  	PRIMARY KEY(bioentry_id,reference_id,rank)
);

CREATE INDEX bioentryref_ref ON bioentry_reference(reference_id);

CREATE TABLE comment (
  	comment_id  	INTEGER PRIMARY KEY,
  	bioentry_id    	INTEGER,
  	comment_text   	TEXT NOT NULL,
  	rank   		SMALLINT NOT NULL DEFAULT 0,
  	UNIQUE(bioentry_id, rank)
);

CREATE TABLE bioentry_qualifier_value (
	bioentry_id   		INTEGER,
   	term_id  		INTEGER,
   	value         		TEXT,
	rank			INT(5) NOT NULL DEFAULT 0,
	UNIQUE (bioentry_id,term_id,rank)
);

CREATE INDEX bioentryqual_trm ON bioentry_qualifier_value(term_id);

CREATE TABLE seqfeature (
   	seqfeature_id 		INTEGER PRIMARY KEY,
   	bioentry_id   		INTEGER,
   	type_term_id		INTEGER,
   	source_term_id  	INTEGER,
	display_name		VARCHAR(64),
   	rank 			SMALLINT  NOT NULL DEFAULT 0,
	UNIQUE (bioentry_id,type_term_id,source_term_id,rank)
);

CREATE INDEX seqfeature_trm  ON seqfeature(type_term_id);
CREATE INDEX seqfeature_fsrc ON seqfeature(source_term_id);

CREATE TABLE seqfeature_relationship (
        seqfeature_relationship_id INTEGER PRIMARY KEY,
   	object_seqfeature_id	INTEGER,
   	subject_seqfeature_id 	INTEGER,
   	term_id 	        INTEGER,
   	rank 			INT(5),
	UNIQUE (object_seqfeature_id,subject_seqfeature_id,term_id)
);

CREATE INDEX seqfeaturerel_trm   ON seqfeature_relationship(term_id);
CREATE INDEX seqfeaturerel_child ON seqfeature_relationship(subject_seqfeature_id);

CREATE TABLE seqfeature_path (
   	object_seqfeature_id	INTEGER,
   	subject_seqfeature_id 	INTEGER,
   	term_id 		INTEGER,
	distance	     	INT(10) ,
	UNIQUE (object_seqfeature_id,subject_seqfeature_id,term_id,distance)
);

CREATE INDEX seqfeaturepath_trm   ON seqfeature_path(term_id);
CREATE INDEX seqfeaturepath_child ON seqfeature_path(subject_seqfeature_id);

CREATE TABLE seqfeature_qualifier_value (
	seqfeature_id 		INTEGER,
   	term_id 		INTEGER,
   	rank 			SMALLINT NOT NULL DEFAULT 0,
   	value  			TEXT NOT NULL,
   	PRIMARY KEY (seqfeature_id,term_id,rank)
);

CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value(term_id);
   
CREATE TABLE seqfeature_dbxref ( 
       	seqfeature_id      INTEGER,
       	dbxref_id          INTEGER,
  	rank  		   SMALLINT,
	PRIMARY KEY (seqfeature_id,dbxref_id)
);

CREATE INDEX feadblink_dbx  ON seqfeature_dbxref(dbxref_id);

CREATE TABLE location (
	location_id		INTEGER PRIMARY KEY,
   	seqfeature_id		INTEGER,
	dbxref_id		INT(10),
	term_id			INT(10),
   	start_pos              	INT(10),
   	end_pos                	INT(10),
   	strand             	TINYINT NOT NULL DEFAULT 0,
   	rank          		SMALLINT NOT NULL DEFAULT 0,
   	UNIQUE (seqfeature_id, rank)
);

CREATE INDEX seqfeatureloc_start ON location(start_pos, end_pos);
CREATE INDEX seqfeatureloc_dbx   ON location(dbxref_id);
CREATE INDEX seqfeatureloc_trm   ON location(term_id);

CREATE TABLE location_qualifier_value (
	location_id		INTEGER,
   	term_id 		INTEGER,
   	value  			VARCHAR(255) NOT NULL,
   	int_value 		INT(10),
	PRIMARY KEY (location_id,term_id)
);

CREATE INDEX locationqual_trm ON location_qualifier_value(term_id);

-- SQLite does not enforce foreign key constraints. There are some trigger
-- based ways to replicate this:
--
-- http://www.sqlite.org/cvstrac/wiki?p=ForeignKeyTriggers
--
-- Currently no foreign key constraints are added.
