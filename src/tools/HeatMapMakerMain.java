package tools;

import algo.FullHeatMap;
import algo.FullHeatMapXML;
import com.sun.org.apache.xml.internal.serialize.XMLSerializer;
import io.IOUtils;
import org.w3c.dom.Document;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class HeatMapMakerMain extends Tool {
    public static final String NAME = "heatmap-maker";
    public static final String DESCRIPTION = "constructs heatmap with dendrogram for distance matrix";


    public final Parameter<File> matrixFile = addParameter(new FileParameterBuilder("matrix-file")
            .mandatory()
            .withShortOpt("i")
            .withDescription("file with distance matrix")
            .create());

    public final Parameter<Boolean> withoutRenumbering = addParameter(new BoolParameterBuilder("without-renumbering")
            .important()
            .withShortOpt("wr")
            .withDescription("don't renumber samples in the heatmap")
            .withDescriptionShort("Don't renumber samples")
            .withDescriptionRu("Не перенумеровывать образцы в итоговой матрице и в тепловой карте")
            .withDescriptionRuShort("Не перенумеровывать образцы")
            .create());

    public final Parameter<File> newMatrixFile = addParameter(new FileParameterBuilder("newMatrix-file")
            .optional()
            .withDefaultComment("<dist-matrix-file>_renumbered.txt")
            .withDescription("resulting renumbered matrix file")
            .create());

    public final Parameter<File> heatmapFile = addParameter(new FileParameterBuilder("heatmap-file")
            .optional()
            .withDefaultComment("<dist-matrix-file>_heatmap.png")
            .withDescription("resulting heatmap file")
            .create());

    public final Parameter<Boolean> invertColors = addParameter(new BoolParameterBuilder("invert-colors")
            .withDescription("invert colors in heatmap")
            .create());

    public final Parameter<String> outputFormat = addParameter(new StringParameterBuilder("output-format")
            .withDefaultValue("%.4f")
            .withDescription("output format for distance values")
            .create());

    public File[] outputDescFiles = null;


    private double[][] matrix;
    private String[] names;


    // output parameters
    private final InMemoryValue<File> heatmapFilePr = new InMemoryValue<File>();
    public final InValue<File> heatmapFileOut = addOutput("heatmap-file", heatmapFilePr, File.class);
    private final InMemoryValue<File> newMatrixFilePr = new InMemoryValue<File>();
    public final InValue<File> newMatrixFileOut = addOutput("newMatrix-file-out", newMatrixFilePr, File.class);



    @Override
    protected void runImpl() throws ExecutionFailedException {
        // parsing input matrix...
        try {
            parseMatrix(matrixFile.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't read matrix file " + matrixFile.get(), e);
        } catch (ParseException e) {
            throw new ExecutionFailedException("Can't read matrix file " + matrixFile.get(), e);
        }


        // creating full heat map
        FullHeatMap maker = new FullHeatMap(matrix, names, invertColors.get());
        BufferedImage image = maker.createFullHeatMap(!withoutRenumbering.get());

        Document document = new FullHeatMapXML(matrix, names, invertColors.get()).createFullHeatMap(!withoutRenumbering.get());

        // saving results
        String filePrefix = FileUtils.removeExtension(matrixFile.get().getPath(), ".txt");

        String newMatrixPath = filePrefix + "_renumbered.txt";
        if (newMatrixFile.get() != null) {
            newMatrixPath = newMatrixFile.get().getPath();
            newMatrixPath = newMatrixPath.replace("$DT", startTimestamp);
        }
        if (withoutRenumbering.get()) {
            // OK
            newMatrixPath = matrixFile.get().getPath();
        } else {
            // should renumber
            try {
                DistanceMatrixCalculatorMain.printMatrix(matrix, newMatrixPath, names, maker.perm, outputFormat.get());
            } catch (FileNotFoundException e) {
                throw new ExecutionFailedException("Can't save renumbered matrix to file " + newMatrixPath, e);
            }
            info("Renumbered matrix saved to " + newMatrixPath);
        }
        newMatrixFilePr.set(new File(newMatrixPath));


        String heatmapPath = FileUtils.removeExtension(newMatrixPath, ".txt") + "_heatmap.png";
        if (heatmapFile.get() != null) {
            heatmapPath = heatmapFile.get().getPath();
            heatmapPath = heatmapPath.replace("$DT", startTimestamp);
        }
        try {
            ImageIO.write(image, "png", new File(heatmapPath));
            XMLSerializer xmlSerializer = new XMLSerializer(new FileWriter(FileUtils.removeExtension(heatmapPath, ".png") + ".svg"), null);
            xmlSerializer.serialize(document);
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't save image to file " + heatmapPath, e);
        }
        info("Heatmap for matrix saved to " + heatmapPath);
        heatmapFilePr.set(new File(heatmapPath));
    }


    private void parseMatrix(File f) throws IOException, ExecutionFailedException, ParseException {
        debug("Parsing matrix from file " + f);

        BufferedReader in = new BufferedReader(new FileReader(f));
        List<String> data = new ArrayList<String>();
        while (in.ready()) {
            data.add(in.readLine());
        }
        in.close();

        if (data.size() == 0) {
            throw new ExecutionFailedException("No data to read in matrix file " + f);
        }
        StringTokenizer st = new StringTokenizer(data.get(0), DistanceMatrixCalculatorMain.SEPARATOR);
        int fn = st.countTokens();
        int sn = data.size();

        if (fn > sn) {
            throw new ExecutionFailedException("Can't parse matrix, columns' number > rows' number");
        }
        // fn <= sn

        // splitting matrix into cells
        String[][] dataArray = new String[fn][fn];
        for (int i = 0; i < fn; i++) {
            st = new StringTokenizer(data.get(i), DistanceMatrixCalculatorMain.SEPARATOR);
            int j = 0;
            while (j < fn && st.hasMoreTokens()) {
                dataArray[i][j] = st.nextToken();
                j++;
            }
            if (j != fn || st.hasMoreTokens()) {
                throw new ExecutionFailedException("Can't parse matrix, columns' number is different for different rows");
            }
        }
        for (int i = fn; i < data.size(); i++) {
            st = new StringTokenizer(data.get(i));
            if (st.hasMoreTokens()) {
                throw new ExecutionFailedException("Can't parse matrix, too much rows");
            }
        }

        // converting to matrix and names
        boolean withNames = false;
        int n = dataArray.length;
        if (dataArray[0][0].equals("#")) {  // with names
            withNames = true;
            n--;
        }
        matrix = new double[n][n];
        names = new String[n];
        for (int i = 0; i < n; i++) {
            if (withNames) {
                names[i] = dataArray[0][i + 1];
            } else {
                names[i] = (i + 1) + " library";
            }
        }
        int dx = withNames ? 1 : 0;
        NumberFormat nf = NumberFormat.getInstance();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = nf.parse(dataArray[i + dx][j + dx]).doubleValue();
            }
        }
        // OK, done
    }


    @Override
    protected void cleanImpl() {
        matrix = null;
        names = null;
    }

    @Override
    protected void postprocessing() {
        if (!withoutRenumbering.get()) {
            IOUtils.tryToAppendDescription(outputDescFiles,
                    newMatrixFileOut.get(),
                    "File with resulted distance matrix between samples with new order based on adjacency of the samples"
            );
        }
        IOUtils.tryToAppendDescription(outputDescFiles,
                heatmapFileOut.get(),
                "Image file with heatmap and dendrogram between samples"
        );
    }


    public HeatMapMakerMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new HeatMapMakerMain().mainImpl(args);
    }
}
