package tools;

import algo.FullHeatMap;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.SimpleDateFormat;
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
            .withDescription("do not renumber samples in the heatmap")
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
        }


        // creating full heat map
        FullHeatMap maker = new FullHeatMap(matrix, names);
        BufferedImage image = maker.createFullHeatMap(!withoutRenumbering.get());



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
                DistanceMatrixCalculatorMain.printMatrix(matrix, newMatrixPath, names, maker.perm);
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
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't save image to file " + heatmapPath, e);
        }
        info("Heatmap for matrix saved to " + heatmapPath);
        heatmapFilePr.set(new File(heatmapPath));
    }


    private void parseMatrix(File f) throws IOException, ExecutionFailedException {
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
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = Double.parseDouble(dataArray[i + dx][j + dx]);
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
        if (outputDescFiles != null) {
            for (File f : outputDescFiles) {
                try {
                    PrintWriter out = new PrintWriter(new FileWriter(f, true));
                    if (!withoutRenumbering.get()) {
                        out.println();
                        out.println(newMatrixFileOut.get());
                        out.println("   File with resulted distance matrix between samples with new order based on adjacency of the samples");
                    }
                    out.println();
                    out.println(heatmapFileOut.get());
                    out.println("   Image file with heatmap and dendrogram between samples");
                    out.println();
                    out.close();
                } catch (IOException e) {
                    // does not matter
                }
            }
        }
    }


    public HeatMapMakerMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new HeatMapMakerMain().mainImpl(args);
    }
}
