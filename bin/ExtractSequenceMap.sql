/* Set the working database */;
use SchlossHanniganSequencing;

/* Get the adapter information */;
SELECT * 
FROM SampleTracking, IndexPrimer001, IndexPrimer002 
WHERE SampleTracking.Index1=IndexPrimer001.Index1 
      AND SampleTracking.Index2=IndexPrimer002.Index2 
      AND SampleTracking.SchlossID="NexteraXT002" 
INTO OUTFILE '/Users/Hannigan/git/Hannigan-2016-ColonCancerVirome/data/NexteraXT002Map.tsv';
