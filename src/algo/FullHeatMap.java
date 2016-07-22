package algo;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * Constructing full heat map with dendrogram for n objects.
 */
public class FullHeatMap {

    // ===========================   Fixed params   ===========================

    final int cellSize = 40;
    final int borderLineSize = 2;
    final int gridSize;

    final int dx_before_dendrogram = 100;
    final int dx_for_dendrogram;
    final int dx_for_names = 350;
    final int dy_for_names = 200;
    final int dy_before_color = 50;
    final int dy_after_color = 100;

    final int dx_scale_text = 120;
    final int dx_scale = 300;
    final int dy_scale = 30;


    public Color lowColor = Color.WHITE;
    public Color highColor = new Color(12, 61, 138);
    public boolean drawInnerLines = true;
    public Color innerLinesColor = Color.GRAY;

    final Font font = new Font(Font.SANS_SERIF, Font.BOLD, 16);


    // ===========================   Variables   ==============================

    public final int n;
    public final String[] names;
    /**
     * Distance matrix with values from 0 to 1
     */
    protected final double[][] distMatrix;
    protected final double low, high;

    public final int[] perm;


    public FullHeatMap(double[][] distMatrix, String[] names) {
        this(distMatrix, 0, 1, names, false);
    }

    public FullHeatMap(double[][] distMatrix, String[] names, boolean invertColors) {
        this(distMatrix, 0, 1, names, invertColors);
    }

    public FullHeatMap(double[][] distMatrix, double low, double high, String[] names, boolean invertColors) {
        n = distMatrix.length;
        this.distMatrix = distMatrix;
        this.low = low;
        this.high = high;
        this.names = names;
        gridSize = n * cellSize + (n + 1) * borderLineSize;
        dx_for_dendrogram = getDendrogramSize(n);
        perm = new int[n];
        for (int i = 0; i < n; i++) {
            perm[i] = i;
        }
        if (invertColors) {
            Color tmp = highColor;
            highColor = lowColor;
            lowColor = tmp;
            innerLinesColor = Color.LIGHT_GRAY;
        }
    }


    static int getDendrogramSize(int n) {
        if (n < 2) {
            return 0;
        }
        int cn = 2;
        int cs = 50;
        while (cn < n) {
            cn++;
            if (3 <= cn && cn <= 5) {
                cs += 50;
            }
            if (6 <= cn && cn <= 10) {
                cs += 30;
            }
            if (11 <= cn && cn <= 20) {
                cs += 15;
            }
        }
        return cs;
    }


    // ==========================   Drawing heat map   ============================

    public BufferedImage createHeatMap() {
        BufferedImage image = new BufferedImage(gridSize, gridSize, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D graphics = image.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.fillRect(0, 0, gridSize, gridSize);

        // drawing grid
        graphics.setColor(Color.GRAY);
        graphics.setStroke(new BasicStroke(borderLineSize));
        int left = borderLineSize / 2;
        int right = n * (cellSize + borderLineSize) + borderLineSize / 2;
        // drawing border
        graphics.drawLine(0, left, gridSize - 1, left);
        graphics.drawLine(left, 0, left, gridSize - 1);
        graphics.drawLine(0, right, gridSize - 1, right);
        graphics.drawLine(right, 0, right, gridSize - 1);
        if (drawInnerLines) {
            graphics.setColor(innerLinesColor);
            for (int i = 1; i < n; i++) {
                int coord = i * (cellSize + borderLineSize) + borderLineSize / 2;
                graphics.drawLine(borderLineSize+1, coord, gridSize - 1 - borderLineSize, coord);
                graphics.drawLine(coord, borderLineSize+1, coord, gridSize - 1 - borderLineSize);
            }
        }

        // filling
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                graphics.setColor(getColor(distMatrix[perm[j]][perm[i]]));
                graphics.fillRect(i * cellSize + (i + 1) * borderLineSize, j * cellSize + (j + 1) * borderLineSize,
                        cellSize, cellSize);
            }
        }

        return image;
    }

    public BufferedImage createHeatMapWithNames() {
        int width = gridSize + dx_for_names; // +dx for row names
        int height = gridSize + dy_for_names;

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D graphics = image.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.fillRect(0, 0, width, height);

        // drawing heat map itself
        BufferedImage heatMap = createHeatMap();
        graphics.drawImage(heatMap, 0, dy_for_names, null);

        // drawing row and column names
        graphics.setColor(Color.BLACK);
        graphics.setFont(font);
        AffineTransform originalTransform = graphics.getTransform();
        AffineTransform newTransform = new AffineTransform(originalTransform);
        newTransform.rotate(-Math.PI / 2, gridSize / 2.0, dy_for_names + gridSize / 2.0);
        for (int i = 0; i < n; i++) {
            int yc = dy_for_names + (i + 1) * (cellSize + borderLineSize) - cellSize / 3;
            graphics.drawString(names[perm[i]], gridSize + 10, yc); // row name

            graphics.setTransform(newTransform);
            graphics.drawString(names[perm[i]], gridSize + 10, yc); // column name
            graphics.setTransform(originalTransform);
        }

        return image;
    }

    protected Color getColor(double value) {
        double p = (value - low) / (high - low);
        int r = Math.min(255,
                lowColor.getRed() + (int) Math.round(p * (highColor.getRed() - lowColor.getRed())));
        int g = Math.min(255,
                lowColor.getGreen() + (int) Math.round(p * (highColor.getGreen() - lowColor.getGreen())));
        int b = Math.min(255,
                lowColor.getBlue() + (int) Math.round(p * (highColor.getBlue() - lowColor.getBlue())));
        return new Color(r, g, b);
    }



    // ============================   Clustering objects   ============================

    class Node {
        int no = -1;    // object number, or -1 - if not a leaf
        int leafs = 1;  // number of leafs

        Node left, right;
        double distance;

        int dy, dx; // dy - from top to bottom, dx - from right to left
    }


    /**
     * Clusters objects, works in O(n^3) time.
     * @return the resulting root node.
     */
    public Node clusterObjects() {
        // preparing
        Node[] nodes = new Node[n];
        for (int i = 0; i < n; i++) {
            Node c = new Node();
            c.no = i;
            c.dx = 0;
            c.dy = (i + 1) * borderLineSize + i * cellSize + cellSize / 2;
            nodes[i] = c;
        }

        double[][] dist = new double[n][n];
        int[] g1 = new int[n];
        int[] g2 = new int[n];

        for (int i = 0; i < n; i++) {
            getGroup(nodes[i], g1, 0);
            dist[i][i] = 0;
            for (int j = i + 1; j < n; j++) {
                getGroup(nodes[j], g2, 0);
                dist[i][j] = dist[j][i] = distanceBetweenGroups(g1, 1, g2, 1);
            }
        }

        // clustering
        int count = n;
        Node root = (n > 0) ? nodes[0] : null;
        while (count > 1) {
//            System.out.println("Count = " + count);

            double minDist = Double.MAX_VALUE;
            int i = -1, j = -1;

            for (int ii = 0; ii < n; ii++) {
                for (int jj = ii + 1; jj < n; jj++) {
                    if (nodes[ii] != null && nodes[jj] != null && dist[ii][jj] < minDist) {
                        minDist = dist[ii][jj];
                        i = ii;
                        j = jj;
                    }
                }
            }

            if (i == -1 || minDist < 0) {
                throw new RuntimeException("Internal error. Wrong minDist index.");
            }

            // merging node i and j
//            System.out.println("Merging " + i + " and " + j + ", min dist = " + minDist);
            root = new Node();
            root.left = nodes[i];
            root.right = nodes[j];
            root.distance = minDist;
            root.leafs = root.left.leafs + root.right.leafs;

            root.dy = (root.left.dy + root.right.dy) / 2;
            root.dx = (int) (minDist / high * dx_for_dendrogram);

            // updating
            nodes[i] = root;
            nodes[j] = null;
            int g1n = getGroup(root, g1, 0);
            for (int ii = 0; ii < n; ii++) {
                dist[ii][j] = dist[j][ii] = -1;
                if (ii != i) {
                    int g2n = getGroup(nodes[ii], g2, 0);
                    dist[ii][i] = dist[i][ii] = distanceBetweenGroups(g1, g1n, g2, g2n);
                }
            }
            count--;
        }

//        System.out.println("Done");
        return root;
    }

    protected double distanceBetweenGroups(int[] g1, int g1n, int[] g2, int g2n) {
        if (g1n == 0 || g2n == 0) {
            return -1;
        }
        double sum = 0;
        for (int i = 0; i < g1n; i++) {
            for (int j = 0; j < g2n; j++) {
                sum += distMatrix[g1[i]][g2[j]];
            }
        }
        return sum / g1n / g2n;
    }

    int getGroup(Node n, int[] e, int first) {
        if (n == null) {
            return 0;
        }
        if (n.no >= 0) {
            e[first] = n.no;
            return 1;
        }
        int c = getGroup(n.left,  e, first);
        c    += getGroup(n.right, e, first + c);
        return c;
    }

    protected void renumber(Node node, int first) {
        if (node.no >= 0) {
            // node is leaf
            perm[first] = node.no;
            node.dy = (first + 1) * borderLineSize + first * cellSize + cellSize / 2;
            return;
        }

        renumber(node.left, first);
        renumber(node.right, first + node.left.leafs);
        node.dy = (node.left.dy + node.right.dy) / 2;
    }


    public BufferedImage createLeftDendrogram(boolean shouldRenumber) {
        int width = dx_before_dendrogram + dx_for_dendrogram;
        int height = gridSize;
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D graphics = image.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.fillRect(0, 0, width, height);

        Node root = clusterObjects();
        if (shouldRenumber) {
            renumber(root, 0);
        }

        graphics.setColor(Color.BLACK);
        graphics.setStroke(new BasicStroke(borderLineSize));
        drawClusterNode(root, graphics, width);

        return image;
    }

    void drawClusterNode(Node n, Graphics2D graphics, int width) {
        if (n == null || n.no >= 0) {
            return;
        }

        graphics.drawLine(width - n.left.dx, n.left.dy, width - n.dx, n.left.dy);
        graphics.drawLine(width - n.right.dx, n.right.dy, width - n.dx, n.right.dy);
        graphics.drawLine(width - n.dx, n.left.dy, width - n.dx, n.right.dy);

        drawClusterNode(n.left, graphics, width);
        drawClusterNode(n.right, graphics, width);

    }



    // ===========================   Drawing full heat map   ================================

    public BufferedImage createColorScale() {
        int dy_before = 20;
        int dy_after = 30;

        int width = dx_scale_text + dx_scale + 1, height = dy_before + dy_scale + dy_after;

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D graphics = image.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.fillRect(0, 0, width, height);

        graphics.setColor(Color.BLACK);
        graphics.setFont(font);
        graphics.drawString("Color", 10, dy_before + dy_scale / 2 + 10);
        int dy_dist = dy_before + dy_scale + 20;
        graphics.drawString("Distance", 10, dy_dist);

//        graphics.drawRect(dx_text, dy_before, dx_scale, dy_scale);
        int cellSize = dx_scale / 6;
        for (int i = 0; i <= 5; i++) {
            double v = low + (high - low) * i / 5.0;
            graphics.setColor(getColor(v));
            graphics.fillRect(dx_scale_text + i * cellSize, dy_before, cellSize, dy_scale);
            graphics.setColor(Color.BLACK);
            graphics.drawString(String.format("%.1f", v), dx_scale_text + i * cellSize + cellSize / 5, dy_dist);
        }

        return image;
    }

    public BufferedImage createFullHeatMap(boolean shouldRenumber) {
        int width = dx_before_dendrogram + dx_for_dendrogram + gridSize + dx_for_names;
        int height = dy_for_names + gridSize + dy_before_color + dy_scale + dy_after_color;
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D graphics = image.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.fillRect(0, 0, width, height);

        if (n >= 2) {
            BufferedImage dendrogram = createLeftDendrogram(shouldRenumber);
            graphics.drawImage(dendrogram, 0, dy_for_names, null);
        }

        BufferedImage heatMap = createHeatMapWithNames();
        graphics.drawImage(heatMap, dx_before_dendrogram + dx_for_dendrogram, 0, null);

        BufferedImage color = createColorScale();
        graphics.drawImage(color, dx_before_dendrogram + dx_for_dendrogram - 120,
                dy_for_names + gridSize + dy_before_color, null);

        return image;
    }


    public static void main(String[] args) throws IOException {
        double[][] matrix = {{  0, 0.9,   1, 0.9, 0.8},
                             {0.9,   0, 0.5, 0.6, 0.1},
                             {  1, 0.5,   0, 0.2, 0.6},
                             {0.9, 0.6, 0.2,   0, 0.7},
                             {0.8, 0.1, 0.6, 0.7, 0}};
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                for (int k = 0; k < matrix.length; k++) {
//                    if (i != j && i != k && j != k) {
                        if (matrix[i][j] + matrix[j][k] < matrix[i][k]) {
                            System.out.println("Wrong for i = " + i + ", j = " + j + ", k = " + k);
                            System.out.println((matrix[i][j] + matrix[j][k]) + " < " + matrix[i][k]);
                            return;
                        }
//                    }
                }
            }
        }
        String[] names = {"Abracadabra library got from some url", "Boldii fish lib2s",
                "Cnot", "Don't know", "mmmm"};

        FullHeatMap hm = new FullHeatMap(matrix, names);
        BufferedImage image = hm.createFullHeatMap(true);

        ImageIO.write(image, "png", new File("test.png"));
    }

}
