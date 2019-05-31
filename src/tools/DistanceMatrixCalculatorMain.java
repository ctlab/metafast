package tools;

import io.IOUtils;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrixCalculatorMain extends Tool {
    public static final String NAME = "dist-matrix-calculator";
    public static final String DESCRIPTION = "Calculates distance matrix using features values";

    public static final String SEPARATOR = "\t";


    public final Parameter<File[]> featuresFiles = addParameter(new FileMVParameterBuilder("features")
            .mandatory()
            .withDescription("features values files (for the same components)")
            .create());

    public final Parameter<Boolean> withoutNames = addParameter(new BoolParameterBuilder("without-names")
            .withShortOpt("wn")
            .optional()
            .withDescription("do not print matrix row and column names (as given file names)")
            .create());

    public final Parameter<File> matrixFile = addParameter(new FileParameterBuilder("matrix-file")
            .optional()
            .withDefaultValue(workDir.append("dist_matrix_$DT_original_order.txt"))
            .withDefaultComment("<workDir>/dist_matrix_<date>_<time>_original_order.txt")
            .withDescription("resulting distance matrix file")
            .create());

    public final Parameter<String> outputFormat = addParameter(new StringParameterBuilder("output-format")
            .withDefaultValue("%.4f")
            .withDescription("output format for distance values")
            .create());

    public File[] outputDescFiles = null;


    @Override
    protected void runImpl() throws ExecutionFailedException {
        List<List<Double>> features = new ArrayList<List<Double>>();

        for (File featuresFile : featuresFiles.get()) {
            try {
                features.add(readVector(featuresFile));
            } catch (IOException e) {
                throw new ExecutionFailedException("Failed to read features from " + featuresFile);
            }
        }

        int cnt = features.size();
        double[][] distMatrix = new double[cnt][cnt];

        for (int i = 0; i < cnt; i++) {
            for (int j = i + 1; j < cnt; j++) {
                distMatrix[i][j] = brayCurtisDistance(features.get(i), features.get(j));
                distMatrix[j][i] = distMatrix[i][j];
            }
        }

        String matrixPath = matrixFile.get().getPath().replace("$DT", startTimestamp);
        String[] names = null;
        if (!withoutNames.get()) {
            names = new String[featuresFiles.get().length];
            for (int i = 0; i < names.length; i++) {
                names[i] = featuresFiles.get()[i].getName();
                names[i] = FileUtils.removeExtension(names[i], "vec");
            }
        }

        try {
            printMatrix(distMatrix, matrixPath, names, null, outputFormat.get());
            info("Distance matrix printed to " + matrixPath);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Failed to print matrix to " + matrixPath);
        }
        matrixFile.set(new File(matrixPath));
    }

    public static void printMatrix(double[][] matrix, String fp, String[] names, int[] perm, String format) throws FileNotFoundException {
        File f = new File(fp);
        FileUtils.makeSubDirsOnly(f);
        PrintWriter out = new PrintWriter(f);

        if (names != null) {
            out.print("#");
            for (int i = 0; i < names.length; i++) {
                out.print(SEPARATOR);
                out.print(names[(perm == null) ? i : perm[i]]);
            }
            out.println();
        }

        for (int i = 0; i < matrix.length; i++) {
            if (names != null) {
                out.print(names[(perm == null) ? i : perm[i]] + SEPARATOR);
            }
            for (int j = 0; j < matrix[i].length; j++) {
                if (j > 0) {
                    out.print(SEPARATOR);
                }
                if (perm == null) {
                    out.printf(format, matrix[i][j]);
                } else {
                    out.printf(format, matrix[perm[i]][perm[j]]);
                }
            }
            out.println();
        }

        out.close();
    }

    private List<Double> readVector(File file) throws IOException {
        List<Double> ans = new ArrayList<Double>();

        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.isEmpty()) {
                ans.add(Double.parseDouble(line));
            }
        }
        reader.close();

        return ans;
    }

    private double brayCurtisDistance(List<Double> vector1, List<Double> vector2) {
        assert vector1.size() == vector2.size();

        double sumdiff = 0, sum = 0;

        for (int pos = 0; pos < vector1.size(); pos++) {
            sumdiff += Math.abs(vector1.get(pos) - vector2.get(pos));
            sum += Math.abs(vector1.get(pos)) + Math.abs(vector2.get(pos));
        }

        assert sum > 0;
        return sumdiff / sum;
    }


    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void postprocessing() {
        IOUtils.tryToAppendDescription(outputDescFiles,
                matrixFile.get(),
                "File with resulted distance matrix between samples keeping original order"
        );
    }


    public static void main(String[] args) {
        new SeqBuilderMain().mainImpl(args);
    }

    public DistanceMatrixCalculatorMain() {
        super(NAME, DESCRIPTION);
    }
}
