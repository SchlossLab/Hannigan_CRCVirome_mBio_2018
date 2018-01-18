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
