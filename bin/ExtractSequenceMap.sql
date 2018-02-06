/* Set the working database */;
use SchlossHanniganSequencing;

/* Get the adapter information */;
SELECT * 
FROM SampleTracking, IndexPrimer001, IndexPrimer002, PureViromeTracking 
WHERE SampleTracking.Index1=IndexPrimer001.Index1
      AND SampleTracking.Index2=IndexPrimer002.Index2
      AND SampleTracking.SchlossID="NexteraXT003"
      AND SampleTracking.UniqueSampleID=PureViromeTracking.UniqueSampleID
INTO OUTFILE '/Users/Hannigan/git/Hannigan-2016-ColonCancerVirome/data/metadata/NexteraXT003Map.tsv';

SELECT * 
FROM SampleTracking, IndexPrimer001, IndexPrimer002, PureViromeTracking 
WHERE SampleTracking.Index1=IndexPrimer001.Index1
      AND SampleTracking.Index2=IndexPrimer002.Index2
      AND SampleTracking.SchlossID="NexteraXT004"
      AND SampleTracking.UniqueSampleID=PureViromeTracking.UniqueSampleID
INTO OUTFILE '/Users/Hannigan/git/Hannigan-2016-ColonCancerVirome/data/metadata/NexteraXT004Map.tsv';

/* Also get the headers to make a human readable table */;
SELECT
	'UniqueSampleID',
	'MetaID',
	'Investigator',
	'LibraryDate',
	'Index1',
	'Index2',
	'Platform',
	'SchlossID',
	'PcrCycles',
	'QubitAfterPrep',
	'TagmentationTime',
	'AvgFragmentLength',
	'CoreSampleNumber',
	'Index1',
	'Barcode1',
	'RevCompBarcode1',
	'Adapter1',
	'Index2',
	'Barcode2',
	'Adapter2',
	'UniqueSampleID',
	'SampleID',
	'ProjectID',
	'StudyID',
	'Processor',
	'QubitConc',
	'PureDate',
	'HostOrg',
	'TimePoint',
	'Disease'
	UNION
SELECT * 
FROM SampleTracking, IndexPrimer001, IndexPrimer002, PureViromeTracking 
WHERE SampleTracking.Index1=IndexPrimer001.Index1
      AND SampleTracking.Index2=IndexPrimer002.Index2
      AND SampleTracking.SchlossID="NexteraXT003"
      AND SampleTracking.UniqueSampleID=PureViromeTracking.UniqueSampleID
INTO OUTFILE '/Users/Hannigan/git/Hannigan-2016-ColonCancerVirome/data/metadata/HumanRead-Virome-NexteraXT003Map.tsv';

SELECT
	'UniqueSampleID',
	'MetaID',
	'Investigator',
	'LibraryDate',
	'Index1',
	'Index2',
	'Platform',
	'SchlossID',
	'PcrCycles',
	'QubitAfterPrep',
	'TagmentationTime',
	'AvgFragmentLength',
	'CoreSampleNumber',
	'Index1',
	'Barcode1',
	'RevCompBarcode1',
	'Adapter1',
	'Index2',
	'Barcode2',
	'Adapter2',
	'UniqueSampleID',
	'SampleID',
	'ProjectID',
	'StudyID',
	'Processor',
	'QubitConc',
	'PureDate',
	'HostOrg',
	'TimePoint',
	'Disease'
	UNION
SELECT * 
FROM SampleTracking, IndexPrimer001, IndexPrimer002, PureViromeTracking 
WHERE SampleTracking.Index1=IndexPrimer001.Index1
      AND SampleTracking.Index2=IndexPrimer002.Index2
      AND SampleTracking.SchlossID="NexteraXT004"
      AND SampleTracking.UniqueSampleID=PureViromeTracking.UniqueSampleID
INTO OUTFILE '/Users/Hannigan/git/Hannigan-2016-ColonCancerVirome/data/metadata/HumanRead-WholeMetagenome-NexteraXT004Map.tsv';
