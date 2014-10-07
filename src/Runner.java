import ru.ifmo.genetics.utils.tool.Tool;

import java.io.*;
import java.lang.reflect.Method;
import java.net.URL;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import static ru.ifmo.genetics.utils.TextUtils.fit;

/**
 * Based on ru.ifmo.genetics.Runner
 */
public class Runner {
    final String TOOLS_PACKAGE = "tools";
    final String SCANNER_PATH = "Runner.class";
    final String USED_FS = "/";

    List<String> classes;
    List<Tool> tools;

    final PrintStream out = System.out;

    void run(String[] args) {
//        args = new String[] {"-t", "antibody-sequences-finder"};
//        System.err.println("Searching for tools...");

        findClasses();
        identifyTools();

//        System.err.println("Found " + tools.size() + " tools");


        boolean printHelp = (args.length == 0) || (args[0].equals("-h") || args[0].equals("--help"));
        if (printHelp) {
            out.println("Fast metagenome analysis toolkit (version 0.1.0)");
            out.println("");
            out.println("Usage: java [<JVM options>] -jar <path-to-jar>/metafast.jar [<Launch options>] [<Tool parameters>]");
            out.println("");
            out.println("This toolkit allows you to run different tools from it.");
            out.println("To see available tools:              java -jar <path-to-jar>/metafast.jar -ts");
            out.println("To see help for selected tool:       java -jar <path-to-jar>/metafast.jar -t <tool-name>");
            out.println("To run selected tool:                java -jar <path-to-jar>/metafast.jar -t <tool-name> <tool-parameters>");
            out.println("");
            return;
        }

        boolean printTools = (args.length > 0) && (args[0].equals("-ts") || args[0].equals("--tools"));
        if (printTools) {
            out.println("Available tools:");
            out.println();
            for (Tool t : tools) {
                out.println(fit(t.name, 30) + " " + fit(t.getClass().getName(), 40) + " " + t.description);
            }
            out.println();
            return;
        }

        boolean toolIsSet = (args.length > 0) && (args[0].equals("-t") || args[0].equals("--tool"));
        if (toolIsSet) {
            if (args.length < 2) {
                throw new RuntimeException("Tool name isn't specified");
            }

            Tool toolInst = null;
            for (Tool tool : tools) {
                if (tool.name.equals(args[1])) {
                    toolInst = tool;
                }
            }
            if (toolInst == null) {
                throw new RuntimeException("Tool name is incorrect");
            }

            args = Arrays.copyOfRange(args, 2, args.length);
            toolInst.mainImpl(args);
            return;
        }

        out.println("Runner: unknown parameters");
    }

    private void findClasses() {
        classes = new ArrayList<String>();

        URL url = Runner.class.getResource(SCANNER_PATH);
        String protocol = url.getProtocol();
        String path = url.getPath();
        path = path.substring(0, path.length() - SCANNER_PATH.length());

        if (protocol.equals("file")) {
            File toolsDir = new File(path + TOOLS_PACKAGE);
            recursiveWalk(toolsDir);
        } else if (protocol.equals("jar")) {
            path = path.substring(5, path.length() - 2);    // removing 'file:' and '!/'
            scanJarFile(path);
        } else {
            throw new RuntimeException("Unknown protocol " + protocol + ", can't scan files");
        }
    }

    private void recursiveWalk(File dir) {
        File[] files = dir.listFiles();
        assert files != null;
        for (File file : files) {
            if (file.isDirectory()) {
                recursiveWalk(file);
            } else {
                String fileName = file.getAbsolutePath();
                fileName = fileName.replace(File.separator, USED_FS);
                if (isGoodClassFile(fileName)) {
                    classes.add(fileName);
                }
            }
        }
    }

    private void scanJarFile(String jarFileName) {
        try {
            ZipFile jarFile = new ZipFile(new File(jarFileName));
            Enumeration<? extends ZipEntry> entries = jarFile.entries();

            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                String fileName = entry.getName();
                if (fileName.startsWith(TOOLS_PACKAGE) && isGoodClassFile(fileName)) {
                    classes.add(fileName);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private boolean isGoodClassFile(String fileName) {
        int li = fileName.lastIndexOf(USED_FS);
        fileName = (li != -1) ? fileName.substring(li + 1) : fileName;
        int pi = fileName.lastIndexOf('.');
        String type = (pi != -1) ? fileName.substring(pi + 1) : "";
        return type.equals("class") && !fileName.contains("$");
    }

    private void identifyTools() {
        tools = new ArrayList<Tool>();

        ClassLoader cl = ClassLoader.getSystemClassLoader();
        for (String classFileName : classes) {
            String className = getClassName(classFileName);
            Class clazz;
            try {
                clazz = cl.loadClass(className);
            } catch (ClassNotFoundException e) {
                throw new RuntimeException(e);
            }

            if (Tool.class.isAssignableFrom(clazz) && (getMainMethod(clazz) != null)) {
                // clazz extends Tool and has main method
                Tool toolClass = null;
                try {
                    toolClass = (Tool) clazz.newInstance();
                } catch (InstantiationException e) {
                    throw new RuntimeException(e);
                } catch (IllegalAccessException e) {
                    throw new RuntimeException(e);
                }

                tools.add(toolClass);
            }
        }
    }

    private Method getMainMethod(Class clazz) {
        Method main = null;
        try {
            main = clazz.getMethod("main", String[].class);
        } catch (NoSuchMethodException e) {
        }
        return main;
    }

    private String getClassName(String fileName) {
        int si = fileName.indexOf(TOOLS_PACKAGE);
        if (si == -1) {
            throw new RuntimeException("si = -1");
        }
        String className = fileName.substring(si);

        int pi = className.lastIndexOf('.');
        className = className.substring(0, pi);     // removing ".class"
        className = className.replaceAll(USED_FS, ".");

        return className;
    }

    public static void main(String[] args) {
        new Runner().run(args);
    }
}
