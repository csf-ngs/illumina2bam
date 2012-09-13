/*
 * Copyright (C) 2011 GRL
 *
 * This library is free software. You can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


package uk.ac.sanger.npg.picard;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.util.FastqQualityFormat;
import net.sf.picard.util.SolexaQualityConverter;
import net.sf.samtools.SAMUtils;
import net.sf.samtools.util.Iso8601Date;
import net.sf.samtools.util.StringUtil;
/**
 * This class is used decode a mutliplexed fastq file.
 * 
 * 
 * 
 * The read group will be changed and re-added in.
 * 
 * @author gq1@sanger.ac.uk
 * 
 */



public class FastQIndexDecoder extends PicardCommandLine {
  
    private final Log log = Log.getInstance(FastQIndexDecoder.class);
    
    private final String programName = "FastQIndexDecoder ";
    
    private final String programDS = "A command-line tool to decode multiplexed fastq file";
   
    @Usage(programVersion= version)
    public final String USAGE = this.getStandardUsagePreamble() + this.programDS + ". ";
   
    @Option(doc="The input fastq file to decode.")
    public File FASTQ1;
    
    @Option( doc="The input fastq file to decode.", optional=true)
    public File FASTQ2;   
    
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file after decoding.", mutex = {"OUTPUT_DIR"} )
    public File OUTPUT;
    
    @Option(doc="The output directory for bam files for each barcode if you want to split the output", mutex = {"OUTPUT"})
    public File OUTPUT_DIR;
    
    @Option(doc="The prefix for bam or sam file when you want to split output by barcodes", mutex = {"OUTPUT"})
    public String OUTPUT_PREFIX;
    
    @Option(doc="The extension name for split file when you want to split output by barcodes: bam or sam", mutex = {"OUTPUT"})
    public String OUTPUT_FORMAT;
    
    @Option(doc="The tag name used to store barcode read in bam records")
    public String BARCODE_TAG_NAME = "BC";

    @Option(doc="Barcode sequence.  These must be unique, and all the same length.", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<String>();

    @Option(doc="Tab-delimited file of barcode sequences, and optionally barcode name and library name.  " +
            "Barcodes must be unique, and all the same length.  Column headers must be 'barcode_sequence', " +
            "'barcode_name', and 'library_name'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Option(doc="Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Option(doc="Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Option(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Option(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Option(shortName="RG", doc="ID used to link RG header record with RG tag in SAM record, default 1.", optional=true)
    public String READ_GROUP_ID = "1";
    
    @Option(shortName="SM", doc="The name of the sequenced sample, using library name if not given.", optional=true)
    public String SAMPLE_ALIAS;

    @Option(shortName="LB", doc="The name of the sequenced library, default unknown.", optional=true)
    public String LIBRARY_NAME = "unknown";

    @Option(shortName="ST", doc="The name of the study.", optional=true)
    public String STUDY_NAME;

    @Option(shortName="PU", doc="The platform unit, using runfolder name plus lane number if not given.", optional=true)
    public String PLATFORM_UNIT;

    @Option(doc="The start date of the run, read from config file if not given.", optional=true)
    public Iso8601Date RUN_START_DATE;

    @Option(shortName="SC", doc="Sequence center name, default SC for Sanger Center.", optional=true)
    public String SEQUENCING_CENTER = "SC";
    
    @Option(doc="The name of the sequencing technology that produced the read, default ILLUMINA.", optional=true)
    public String PLATFORM = "ILLUMINA";
    
    @Option(shortName="V", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 " +
            "style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled " +
            "scores with a character shift of 33.")
    public FastqQualityFormat QUALITY_FORMAT;
    
    private int barcodeLength;
    
    private IndexDecoder indexDecoder;
    
    private SAMFileWriter out;
    private HashMap<String, SAMFileWriter> outputList;
    private HashMap<String, String> barcodeNameList;
    
    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();
    
    public FastQIndexDecoder() {
    }

    @Override
    protected int doWork() {
        
        String readGroup = READ_GROUP_ID;
        
        this.log.info("Checking input and output file");
        IoUtil.assertFileIsReadable(FASTQ1);
        if(FASTQ2 != null){
            IoUtil.assertFileIsReadable(FASTQ2);
        }        
        if(OUTPUT != null){
            IoUtil.assertFileIsWritable(OUTPUT);
        }
        if(OUTPUT_DIR != null){
            IoUtil.assertDirectoryIsWritable(OUTPUT_DIR);
        }
        IoUtil.assertFileIsWritable(METRICS_FILE);
        
        log.info("Open input file: " + FASTQ1.getName());
        final FastqReader in1  = new FastqReader(FASTQ1);  
      

        final SAMFileHeader header = new SAMFileHeader();
        SAMReadGroupRecord rgr = generateSamReadGroupRecord();
        List<SAMReadGroupRecord> rgrl = new ArrayList<SAMReadGroupRecord>();
        rgrl.add(rgr);
        header.setReadGroups(rgrl);
        
        this.generateOutputFile(header);
                
        log.info("Decoding records");        
        Iterator<FastqRecord> iter1 = in1.iterator();
     
        FastqReader in2 = null;
        Iterator<FastqRecord> iter2 = null;
        boolean isPaired = false;
        if(FASTQ2 != null){
           in2 = new FastqReader(FASTQ2);
           isPaired = true;
        }
        while(iter1.hasNext()){
            
            String barcodeRead = null;

            FastqRecord record = in1.next();            
            String readName = record.getReadHeader();

            Object barcodeReadObject = record.getReadString().substring(0,barcodeLength);
            
            FastqRecord pairedRecord = null;
            if(isPaired){
                pairedRecord = in2.next();
                String readName2 = pairedRecord.getReadHeader();
        // could check if names are identical too complicated now!       
        //        if( !readName.equals(readName2) || !isPaired2 ){
        //            throw new RuntimeException("The paired reads are not together: " + readName + " " + readName2);
         //       }

                Object barcodeReadObject2= pairedRecord.getReadString().substring(0,barcodeLength);
          //      if(barcodeReadObject != null
         //               && barcodeReadObject2 != null
         //               && ! barcodeReadObject.equals(barcodeReadObject2) ){
         //           
         //           throw new RuntimeException("barcode read bases are different in paired two reads: "
         //                   + barcodeReadObject + " " + barcodeReadObject2);
         //       } else if( barcodeRead == null && barcodeReadObject2 != null ){
         //           barcodeRead = barcodeReadObject2.toString();
         //       }                
            }
            barcodeRead = record.getReadString().substring(0,barcodeLength);
            if(barcodeRead == null ){
                throw new RuntimeException("No barcode read found for record: " + readName );
            }
            
            if(barcodeRead.length() < this.barcodeLength){
                throw new RuntimeException("The barcode read length is less than barcode lenght: " + readName );
            }else{            
                barcodeRead = barcodeRead.substring(0, this.barcodeLength);
            }

            IndexDecoder.BarcodeMatch match = this.indexDecoder.extractBarcode(barcodeRead, true); //no quality filtering all are passed
            String barcode = match.barcode;
            
            if( match.matched ) {
               barcode = barcode.toUpperCase();
            } else {
               barcode = "";
            }
            
            SAMRecord samRec1 = new SAMRecord(header);
            SAMRecord samRec2 = null;
  
            String barcodeName = this.barcodeNameList.get(barcode);
            samRec1.setReadUnmappedFlag(true);
            if(match.matched){
              String cutNuc = record.getReadString().substring(0,barcodeLength);
              String cutQual = convertQualityTag(record.getBaseQualityString().substring(0, barcodeLength));
              samRec1.setReadString(record.getReadString().substring(barcodeLength));
              setQuality(samRec1,record.getBaseQualityString().substring(barcodeLength)); 
              samRec1.setAttribute(BARCODE_TAG_NAME, cutNuc);
              samRec1.setAttribute("QT", cutQual);
            }else{
              samRec1.setReadString(record.getReadString());
              setQuality(samRec1,record.getBaseQualityString());  
            }
            samRec1.setReadName(readName + "#" + barcodeName);
            samRec1.setAttribute("RG", readGroup + "#" + barcodeName);
            if (isPaired) {
                samRec2 = new SAMRecord(header);
                samRec2.setReadName(readName + "#" + barcodeName);
                samRec2.setAttribute("RG", readGroup + "#" + barcodeName);
                String cutNuc = "";
                String cutQual = "";
                if(match.matched){
                   cutNuc = pairedRecord.getReadString().substring(0,barcodeLength);
                   cutQual = convertQualityTag(pairedRecord.getBaseQualityString().substring(0, barcodeLength));
                   samRec2.setReadString(pairedRecord.getReadString().substring(barcodeLength));
                   setQuality(samRec2, pairedRecord.getBaseQualityString().substring(barcodeLength));
                }else{
                   samRec2.setReadString(pairedRecord.getReadString());
                   setQuality(samRec2, pairedRecord.getBaseQualityString()); 
                }
                samRec1.setReadPairedFlag(true);
                samRec1.setFirstOfPairFlag(true);
                samRec1.setMateUnmappedFlag(true);
                samRec2.setAttribute(BARCODE_TAG_NAME, cutNuc);
                samRec2.setAttribute("QT", cutQual);
                samRec2.setFirstOfPairFlag(false);
                samRec2.setReadPairedFlag(true);
                samRec2.setSecondOfPairFlag(true);
                samRec2.setMateUnmappedFlag(true);
                samRec2.setReadUnmappedFlag(true);
            }
            
            if( OUTPUT != null ){
                out.addAlignment(samRec1);
                if(isPaired){
                    out.addAlignment(samRec2);
                }
            } else {
                
                SAMFileWriter outPerBarcode = this.outputList.get(barcode);
                outPerBarcode.addAlignment(samRec1);
                if(isPaired){
                    outPerBarcode.addAlignment(samRec2);
                }                
            }
            
        }
        
        if(out != null){
           out.close();
        }
        this.closeOutputList();
        
        log.info("Decoding finished");
        
        
        log.info("Writing out metrhics file");        
        final MetricsFile<IndexDecoder.BarcodeMetric, Integer> metrics = getMetricsFile();        
        indexDecoder.writeMetrics(metrics, METRICS_FILE);
        
        log.info("All finished");

        return 0;
    }
    
    /**
     * Generate read group record
     * 
     * @param platformUnitConfig default platform unit from configure XML, which will be used if not given from command line, and could be null
     * @param runDateConfig default run date from configure XML, which will be used if not given from command line, and could be null
     * @return read group record for BAM header
     */
    public SAMReadGroupRecord generateSamReadGroupRecord(){
        
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord(this.READ_GROUP_ID);
        
        readGroup.setLibrary(this.LIBRARY_NAME);
        
        if(this.SAMPLE_ALIAS == null){
           readGroup.setSample(this.LIBRARY_NAME);
        }else{
            readGroup.setSample(this.SAMPLE_ALIAS);
        }
        
        if( this.STUDY_NAME != null ){
            readGroup.setDescription("Study " + this.STUDY_NAME);
        }
        
        if(this.PLATFORM_UNIT != null){
            readGroup.setPlatformUnit(this.PLATFORM_UNIT);
        }
        
        if(this.RUN_START_DATE != null){
            final Iso8601Date date = new Iso8601Date(RUN_START_DATE);
            readGroup.setRunDate(date);
        }
        
        readGroup.setPlatform(this.PLATFORM);
        
        readGroup.setSequencingCenter(this.SEQUENCING_CENTER);
        
        return readGroup;
    }
    
    public void generateOutputFile(SAMFileHeader header) {
        
        List<IndexDecoder.NamedBarcode> barcodeList = indexDecoder.getNamedBarcodes(); 
        
        this.barcodeNameList = new HashMap<String, String>();
        
        List<SAMReadGroupRecord> oldReadGroupList = header.getReadGroups();        
        List<SAMReadGroupRecord> fullReadGroupList = new ArrayList<SAMReadGroupRecord>();
        
        if (OUTPUT_DIR != null) {
            log.info("Open a list of output bam/sam file per barcode");
            outputList = new HashMap<String, SAMFileWriter>();
        }

        for (int count = 0; count <= barcodeList.size(); count++) {

            String barcodeName = null;
            String barcode = null;
            IndexDecoder.NamedBarcode namedBarcode = null;
            List<SAMReadGroupRecord> readGroupList = new ArrayList<SAMReadGroupRecord>();

            if ( count != 0 ) {
                namedBarcode = barcodeList.get(count - 1);
                barcodeName = namedBarcode.barcodeName;
                barcode = namedBarcode.barcode;
                barcode = barcode.toUpperCase();
            }else{
                barcode = "";
            }

            if (barcodeName == null || barcodeName.equals("")) {
                barcodeName = Integer.toString(count);
            }

            for(SAMReadGroupRecord r : oldReadGroupList){
                    SAMReadGroupRecord newReadGroupRecord = new SAMReadGroupRecord(r.getId() + "#" + barcodeName, r);
                    String pu = newReadGroupRecord.getPlatformUnit();
                    if(pu != null){
                        newReadGroupRecord.setPlatformUnit(pu + "#" + barcodeName);
                    }
                    if(namedBarcode != null){
                        if( namedBarcode.libraryName != null && !namedBarcode.libraryName.equals("") ){
                           newReadGroupRecord.setLibrary(namedBarcode.libraryName);
                        }
                        if( namedBarcode.sampleName !=null && !namedBarcode.sampleName.equals("") ){
                           newReadGroupRecord.setSample(namedBarcode.sampleName);
                        }
                        if(namedBarcode.description != null && !namedBarcode.description.equals("") ){
                            newReadGroupRecord.setDescription(namedBarcode.description);
                        }
                    }
                    readGroupList.add(newReadGroupRecord);
            }
            fullReadGroupList.addAll(readGroupList);


            if (OUTPUT_DIR != null) {
                String barcodeBamOutputName = OUTPUT_DIR
                        + File.separator
                        + OUTPUT_PREFIX
                        + "#"
                        + barcodeName
                        + "."
                        + OUTPUT_FORMAT;
                final SAMFileHeader outputHeader = header.clone();
                outputHeader.setReadGroups(readGroupList);
                this.addProgramRecordToHead(outputHeader, this.getThisProgramRecord(programName, programDS));
                final SAMFileWriter outPerBarcode = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader, true, new File(barcodeBamOutputName));
                outputList.put(barcode, outPerBarcode);
            }
            barcodeNameList.put(barcode, barcodeName);
        }
        
        if (OUTPUT != null) {
            log.info("Open output file with header: " + OUTPUT.getName());
            final SAMFileHeader outputHeader = header.clone();
            outputHeader.setReadGroups(fullReadGroupList);
            this.addProgramRecordToHead(outputHeader, this.getThisProgramRecord(programName, programDS));
            this.out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader, true, OUTPUT);
        }

    }
    
    public void closeOutputList(){
        if( this.outputList != null ){
            for(SAMFileWriter writer: this.outputList.values()){
                writer.close();
            }
        }
    }
    
     /** Based on the type of quality scores coming in, converts them to a numeric byte[] in prhred scale. */
    void convertQuality(final byte[] quals, final FastqQualityFormat version) {
        switch (version)  {
            case Standard:
                SAMUtils.fastqToPhred(quals);
                break ;
            case Solexa:
                solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
                break ;
            case Illumina:
                solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
                break ;
            }
    }

    
    String convertQualityTag(String qualities){
        final byte[] quals = StringUtil.stringToBytes(qualities);
        convertQuality(quals, QUALITY_FORMAT);
        return SAMUtils.phredToFastq(quals);
    }
        
 
    void setQuality(SAMRecord record, String qualities){
        final byte[] quals = StringUtil.stringToBytes(qualities);
        convertQuality(quals, QUALITY_FORMAT);
        record.setBaseQualities(quals);
    }
    
    /**
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        
        final ArrayList<String> messages = new ArrayList<String>();

        if (BARCODE_FILE != null) {
            this.indexDecoder = new IndexDecoder(BARCODE_FILE);
        } else {
            this.indexDecoder = new IndexDecoder(BARCODE);
        }
        
        indexDecoder.setMaxMismatches(this.MAX_MISMATCHES);
        indexDecoder.setMaxNoCalls(MAX_NO_CALLS);
        indexDecoder.setMinMismatchDelta(this.MIN_MISMATCH_DELTA);
        
        indexDecoder.prepareDecode(messages);
        this.barcodeLength = indexDecoder.getBarcodeLength();

        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    public static void main(final String[] argv) {
        System.exit(new FastQIndexDecoder().instanceMain(argv));
    }

}
