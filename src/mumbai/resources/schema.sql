-- This version of schema lacks bin R*tree and supports only queries where
-- the VALUES of bins is supplied as a query parameter. It also has only
-- the indexes revealed to be used by EXPLAIN QUERY PLAN ...

-- Note: ref.name requires bam or other source as this isn't in bai
CREATE TABLE ref (
  id ${ref_id_type} PRIMARY KEY,
  name TEXT NOT NULL UNIQUE
);

-- Note: id is 1-based because it's auto-incrementing integer primary key
CREATE TABLE bam (
  id ${bam_id_type} PRIMARY KEY,
  uuid TEXT NOT NULL UNIQUE,
  filepath TEXT
);

CREATE TABLE chunk (
  bam_id INTEGER NOT NULL,
  ref_id INTEGER NOT NULL,
  bin_id INTEGER NOT NULL,
  chunk_beg ${chunk_chunk_beg_type} NOT NULL,
  chunk_end ${chunk_chunk_end_type} NOT NULL,
  FOREIGN KEY(bam_id) REFERENCES bam(id),
  FOREIGN KEY(ref_id) REFERENCES ref(id)
);

-- linear index
-- Note: that window is zero-based to conform to SAM spec
CREATE TABLE intv (
  bam_id INTEGER NOT NULL,
  ref_id INTEGER NOT NULL,
  -- WINDOW is a keyword in postgres and needs quotes as column name
  "window" INTEGER NOT NULL,
  ioffset ${intv_ioffset_type} NOT NULL,
  FOREIGN KEY(bam_id) REFERENCES bam(id),
  FOREIGN KEY(ref_id) REFERENCES ref(id)
);

-- indexes -- tested good!
CREATE INDEX chunk_bin_ref_bam_idx ON chunk(bin_id, ref_id, bam_id);
CREATE INDEX intv_window_ref_bam_idx ON intv("window", ref_id, bam_id);
