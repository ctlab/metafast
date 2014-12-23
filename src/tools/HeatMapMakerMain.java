package tools;

import algo.FullHeatMap;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
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


    private double[][] matrix;
    private String[] names;


    @Override
    protected void runImpl() throws ExecutionFailedException {
        try {
            parseMatrix(matrixFile.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't read matrix file " + matrixFile.get(), e);
        }

        String fileName = FileUtils.removeExtension(matrixFile.get().toString(), ".txt");
        fileName += "_heatmap.png";

        FullHeatMap maker = new FullHeatMap(matrix, names);
        BufferedImage image = maker.createFullHeatMap();
        try {
            ImageIO.write(image, "png", new File(fileName));
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't save image to file " + fileName, e);
        }

        info("Heatmap for matrix saved to " + fileName);
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

    public HeatMapMakerMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new HeatMapMakerMain().mainImpl(args);
    }
}
