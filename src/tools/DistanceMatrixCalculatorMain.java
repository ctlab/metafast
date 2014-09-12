package tools;

import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ulyantsev on 29.05.14.
 *
 */
public class DistanceMatrixCalculatorMain extends Tool {
    public static final String NAME = "calc-dist";
    public static final String DESCRIPTION = "Calculates distance matrix using features files";


    public final Parameter<File[]> featuresFiles = addParameter(new FileMVParameterBuilder("features")
            .mandatory()
            .withDescription("features files provided from features-builder (for same components)")
            .create());

    public final Parameter<String> separator = addParameter(new StringParameterBuilder("separator")
            .withShortOpt("s")
            .withDescription("matrix values separator")
            .withDefaultValue(" ")
            .withDefaultComment("one white space")
            .create());

    public final Parameter<Boolean> printNames = addParameter(new BoolParameterBuilder("print-names")
            .withShortOpt("pn")
            .optional()
            .withDescription("add to matrix row and column with given file names")
            .withDefaultValue(false)
            .create());

    @Override
    protected void runImpl() throws ExecutionFailedException {
        List<List<Long>> features = new ArrayList<List<Long>>();

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

        String matrixPath = workDir + File.separator + "dist_matrix";
        try {
            printMatrix(distMatrix, matrixPath);
            info("Distance matrix printed to " + matrixPath);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Failed to print matrix to " + matrixPath);
        }
    }

    private void printMatrix(double[][] matrix, String fp) throws FileNotFoundException {
        PrintWriter matrixPW = new PrintWriter(fp);

        if (printNames.get()) {
            matrixPW.print("#");
            for (File featuresFile : featuresFiles.get()) {
                matrixPW.print(separator.get());
                matrixPW.print(featuresFile.getName());
            }
            matrixPW.println();
        }

        for (int i = 0; i < matrix.length; i++) {
            if (printNames.get()) {
                matrixPW.print(featuresFiles.get()[i].getName() + separator.get());
            }
            for (int j = 0; j < matrix[i].length; j++) {
                if (j > 0) {
                    matrixPW.print(separator.get());
                }
                matrixPW.print(matrix[i][j]);
            }
            matrixPW.println();
        }

        matrixPW.close();
    }

    private List<Long> readVector(File file) throws IOException {
        List<Long> ans = new ArrayList<Long>();

        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.isEmpty()) {
                ans.add(Long.parseLong(line));
            }
        }
        reader.close();

        return ans;
    }

    private double brayCurtisDistance(List<Long> vector1, List<Long> vector2) {
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

    public static void main(String[] args) {
        new SeqBuilderMain().mainImpl(args);
    }

    public DistanceMatrixCalculatorMain() {
        super(NAME, DESCRIPTION);
    }
}
