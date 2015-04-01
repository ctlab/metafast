package tools;

import structures.ConnectedComponent;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class CompareReadsAndComponentsMain extends Tool {
    public static final String NAME = "comparison-script";
    public static final String DESCRIPTION = "Statistics";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .mandatory()
            .withShortOpt("cm")
            .withDescription("file with connected components in binary format")
            .create());

    public final Parameter<File> referenceFile = addParameter(new FileParameterBuilder("reference-file")
            .mandatory()
            .withShortOpt("r")
            .withDescription("FASTA file with reference")
            .create());

    public final Parameter<File> samtoolsOutput = addParameter(new FileParameterBuilder("samtools-file")
            .mandatory()
            .withShortOpt("so")
            .withDescription("SamTools view reads output from BAM file")
            .create());

    private List<String> contigsID = new ArrayList<String>();
    private List<String> referenceContigs = new ArrayList<String>();

    private List<int[]> readsBegins = new ArrayList<int[]>();
    private List<int[]> readsEnds = new ArrayList<int[]>();

    private List<ConnectedComponent> components;

    @Override
    protected void runImpl() throws ExecutionFailedException {
        debug("Lets load components");
        components = ConnectedComponent.loadComponents(componentsFile.get());
        debug("Components loaded, loading reference");

        try {
            readReferenceContigs();
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load reference", e);
        }
        debug("Reference loaded");

        try {
            BufferedReader br = new BufferedReader(new FileReader(samtoolsOutput.get()));
            String line;
            int nLines = 0;
            while ((line = br.readLine()) != null) {
                nLines++;
                String[] splitted = line.split("\\s+");

                String id = splitted[2];
                int pos = Integer.parseInt(splitted[3]);
                int readLen = Integer.parseInt(splitted[5].substring(0, splitted[5].length() - 1));
                if (nLines % 1000000 == 0) {
                    debug("lines cnt = " + nLines);
                }

                int idPos = contigsID.indexOf(id);
                if (idPos > -1) {
                    readsBegins.get(idPos)[pos]++;
                    readsEnds.get(idPos)[Math.min(pos + readLen - 1, readsEnds.get(idPos).length - 1)]++;
                }

            }
            br.close();
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load samtools", e);
        }

        debug("Building kmer-to-component map");
        ArrayLong2IntHashMap kmerToComponent =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4);

        for (int compNum = 0; compNum < components.size(); compNum++) {
            for (long kmer : components.get(compNum).kmers) {
                assert kmerToComponent.get(kmer) == 0;
                kmerToComponent.add(kmer, compNum + 1);
            }
        }

        debug("Printing statistics");
        try {
            int[] componentToCount = new int[components.size() + 1];
            PrintWriter pw = new PrintWriter(workDir + File.separator + "reference-to-component");

            long inComponents = 0, inReads = 0, inComponentsAndReads = 0;

            for (int contigNum = 0; contigNum < contigsID.size(); contigNum++) {
                pw.println(contigsID.get(contigNum));

                String contig = referenceContigs.get(contigNum);
                int[] begins = readsBegins.get(contigNum), ends = readsEnds.get(contigNum);
                ShortKmer currentKmer = new ShortKmer(contig.substring(0, k.get()));

                int currentReadsCount = 0;
                for (int pos = 0; pos <= contig.length(); pos++) {
                    if (pos > 0) {
                        currentReadsCount += begins[pos - 1] - (pos > 1 ? ends[pos - 2] : 0);
                    }

                    if (pos >= k.get()) {
                        int componentNum = kmerToComponent.get(currentKmer.toLong());
                        pw.println((pos - k.get()) + " " + componentNum + " " + currentReadsCount);
                        componentToCount[componentNum]++;

                        if (componentNum > 0) {
                            if (currentReadsCount > 0) {
                                inComponentsAndReads++;
                            } else {
                                inComponents++;
                            }
                        } else if (currentReadsCount > 0) {
                            inReads++;
                        }
                    }

                    if (pos < contig.length()) {
                        currentKmer.shiftRight(DnaTools.fromChar(contig.charAt(pos)));
                    }
                }
            }
            pw.close();

            info("just in reads = " + inReads);
            info("just in components = " + inComponents);
            info("in components and reads = " + inComponentsAndReads);

            PrintWriter componentStatPW = new PrintWriter(workDir + File.separator + "components-stat");
            for (int i = 0; i < componentToCount.length; i++) {
                if (componentToCount[i] > 0) {
                    componentStatPW.println(i + " " + componentToCount[i]);
                }
            }
            componentStatPW.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        debug("done");
    }

    protected void readReferenceContigs() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(referenceFile.get()));

        StringBuilder sb = null;
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith(">")) {
                if (sb != null) {
                    referenceContigs.add(sb.toString());
                }
                sb = new StringBuilder();
                contigsID.add(line.substring(1));
            } else {
                assert sb != null;
                sb.append(line);
            }
        }
        if (sb != null) {
            referenceContigs.add(sb.toString().substring(1));
        }

        for (int i = 0; i < contigsID.size(); i++) {
            readsBegins.add(new int[referenceContigs.get(i).length()]);
            readsEnds.add(new int[referenceContigs.get(i).length()]);
        }

        br.close();
    }

    @Override
    protected void cleanImpl() {

    }

    public static void main(String[] args) {
        new CompareReadsAndComponentsMain().mainImpl(args);
    }

    public CompareReadsAndComponentsMain() {
        super(NAME, DESCRIPTION);
    }
}
