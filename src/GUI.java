import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.TextUtils;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import tools.DistanceMatrixBuilderMain;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static ru.ifmo.genetics.utils.TextUtils.NL;

public class GUI extends JPanel {

    DistanceMatrixBuilderMain metafast = new DistanceMatrixBuilderMain();

    // MetaFast params
    Parameter<File[]> inputFiles = metafast.inputFiles;
    JTable inFilesTable;
    JButton addFileButton;
    JButton changeFileButton;
    JButton removeFileButton;


    Parameter<Integer> kmerSize = metafast.k;
    JSpinner kmerSpinner;
    JLabel kmerLabel;

    Parameter<Integer> maxBadFrequency = metafast.maximalBadFrequency;
    JSpinner maxBadFrequencySpinner;
    JLabel maxBadFrequencyLabel;

    Parameter<Integer> minSequenceLength = metafast.minSequenceLength;
    JSpinner minSeqLenSpinner;
    JLabel minSeqLenLabel;

    Parameter<Integer> minComponentSize = metafast.compCutter.minComponentSize;
    JSpinner minCompSizeSpinner;
    JLabel minCompSizeLabel;

    Parameter<Integer> maxComponentSize = metafast.compCutter.maxComponentSize;
    JSpinner maxCompSizeSpinner;
    JLabel maxCompSizeLabel;

    Parameter<Boolean> withoutRenumbering = metafast.heatMapMaker.withoutRenumbering;
    JCheckBox wrCheckBox;


    // Global options
    Parameter<String> workDir = Tool.workDirParameter;
    JTextField workDirField;
    JLabel workDirLabel;
    JButton setWorkDirButton;

    String currentMemory;
    String memory;
    JTextField memoryField;
    JLabel memoryLabel;
    JButton setMemoryButton;

    Parameter<Integer> availableProcessors = Tool.availableProcessorsParameter;
    JSpinner avProcSpinner;
    JLabel avProcLabel;

    boolean currentEA;
    boolean ea;
    JCheckBox eaCheckBox;

    Parameter<Boolean> verbose = Tool.verboseParameter;
    JCheckBox verboseCheckBox;


    Parameter<Boolean> startNow = Parameter.createParameter(new BoolParameterBuilder("start-now")
            .withDefaultValue(false)
            .withDescription("start the analysis now!")
            .create());


    // determine default values
    {
//        inputFiles.set(new File[0]);
        determineMemoryDefaultValue();
        determineEADefaultValue();
    }
    Parameter[] almostAllParams = {workDir, availableProcessors, verbose,
            kmerSize, maxBadFrequency, minSequenceLength,
            minComponentSize, maxComponentSize, withoutRenumbering,
            startNow};
    Parameter[] allParams = {workDir, availableProcessors, verbose,
            kmerSize, maxBadFrequency, minSequenceLength,
            minComponentSize, maxComponentSize, withoutRenumbering,
            startNow, inputFiles};


    // others UI elements
    TitledBorder br1, br2, br3, br4;
    JButton runButton, stopButton;
    int more1State, more2State;
    JButton more1Button, more2Button;

    // status elements
    JLabel runningStageText, runningStage;
    JLabel elapsedTimeText, elapsedTime;
    JLabel remainingTimeText, remainingTime;
    JLabel overallProgressText, overallProgress;
    JLabel logText;
    JTextArea log;


    // some ui text
    int curLang = 0;    // 0=english, 1=russian
    String[][] UIText = {
            /* 0 */      {"1. Choose input files", "1. Выберите файлы с исходными чтениями"},
            /* 1 */      {"2. Set options and parameters", "2. Установите опции и параметры анализа"},
            /* 2 */      {"3. Run the analysis", "3. Запустите процесс анализа"},
            /* 3 */      {"4. Running status", "4. Состояние выполнения"},
            /* 4 */      {"Add files", "Добавить файлы"},
            /* 5 */      {"Change file path", "Выбрать другой файл"},
            /* 6 */      {"Remove selected files", "Удалить выбранные файлы"},
            /* 7 */      {"More", "Больше"},
            /* 8 */      {"Less", "Меньше"},
            /* 9 */      {"Choose", "Выбрать"},
           /* 10 */      {"Set", "Установить"},
           /* 11 */      {"Start the analysis", "Начать процесс анализа"},
           /* 12 */      {"Stop the analysis", "Остановить процесс анализа"},
           /* 13 */      {"No input files were selected.", "Файлы не выбраны."},
           /* 14 */      {"Can't start Java with given memory size:", "Не получается запустить Java-машину с заданным объемом памяти:"},
           /* 15 */      {"Can't set memory while analysis process is running." + NL + "Please, stop the analysis process first.",
                          "Нельзя изменить объем используемой памяти во время выполнения анализа." + NL + "Пожалуйста, сначала остановите процесс анализа, а потом поставьте нужный объем памяти."},
           /* 16 */      {"The analysis is already running!", "Процесс анализа уже идет!"},
           /* 17 */      {"Working directory (WD) contains files from previous run. " + NL + "MetaFast can continue the previous run, " +
                            "or start from the beginning and rewrite old files. " + NL + "What option do you prefer?",
                          "Рабочая директория (WD) содержит файлы от предыдущего запуска. " + NL + "MetaFast может продолжить выполнение предыдущего запуска, " +
                            "или начать анализа с начала и перезаписать старые файлы. " + NL + "Какой вариант использовать?"},
           /* 18 */      {"Start from the beginning and rewrite old files", "Начать анализ с начала и перезаписать старые файлы"},
           /* 19 */      {"Continue the previous run", "Продолжить выполнение прошлого запуска"},
           /* 20 */      {"Do nothing", "Ничего не делать"},
           /* 21 */      {"Analysis process isn't running now!", "Процесс анализа сейчас не выполняется!"},
           /* 22 */      {"No input files were selected!", "Файлы с исходными чтениями не выбраны!"},
           /* 23 */      {"An exception occurred during analysis process:", "Произошла ошибка во время анализа:"},
           /* 24 */      {"The analysis has finished successfully!" + NL + "See working directory and/or output_description.txt file for view the results.",
                          "Процесс анализа завершился успешно!" + NL + "Для просмотра результатов анализа используйте рабочую директорию и/или файл output_description.txt."},
           /* 25 */      {"Please select exactly one file!", "Пожалуйста выберите ровно один файл!"},
           /* 26 */      {"Please select files to remove!", "Пожалуйста выберите файлы для удаления!"},
    };
    
    String[][] UIOptionsText = {
            /* 0 */      {"Memory to use", "Использовать памяти"},
            /* 1 */      {"Enable assertions", "Включить проверку assert'ов"},
    };
    String[][] UIStatusText = {
            /* 0 */      {"Running stage:", "Выполняемый этап:"},
            /* 1 */      {"Elapsed time:", "Затрачено времени:"},
            /* 2 */      {"Remaining time:", "Осталось:"},
            /* 3 */      {"Overall progress:", "Прогресс сборки:"},
            /* 4 */      {"Log", "Лог"},
    };



    public GUI(String[] args) {
        // parsing args
        Tool.parseArgs(Arrays.asList(allParams), args);

        // creating UI
        setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));

        JPanel languagePanel = createLanguagePanel();
        JPanel inputFilesPanel = createInputFilesPanel();
        JPanel optionsPanel = createParamsPanel();
        JPanel controlPanel = createControlPanel();
        JPanel statusPanel = createStatusPanel();

        add(Box.createRigidArea(new Dimension(0, 10)));
        add(languagePanel);
//        add(Box.createRigidArea(new Dimension(0, 10)));
        add(inputFilesPanel);
//        add(Box.createRigidArea(new Dimension(0, 10)));
        add(optionsPanel);
//        add(Box.createRigidArea(new Dimension(0, 10)));
        add(controlPanel);
//        add(Box.createRigidArea(new Dimension(0, 10)));
        add(statusPanel);

        updateUIText();
        updateUIValuesFromGlobalValues();
    }


    JPanel createLanguagePanel() {
        JPanel panel = new JPanel();

        JLabel enLabel = new JLabel("English");
        enLabel.setForeground(Color.BLUE);
        enLabel.setCursor(new Cursor(Cursor.HAND_CURSOR));
        enLabel.setIcon(new ImageIcon(ClassLoader.getSystemResource("images/en.png")));
        enLabel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                curLang = 0;
                updateUIText();
            }
        });

        JLabel ruLabel = new JLabel("Русский");
        ruLabel.setForeground(Color.BLUE);
        ruLabel.setCursor(new Cursor(Cursor.HAND_CURSOR));
        ruLabel.setIcon(new ImageIcon(ClassLoader.getSystemResource("images/ru.png")));
        ruLabel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                curLang = 1;
                updateUIText();
            }
        });

        panel.add(enLabel);
        panel.add(Box.createRigidArea(new Dimension(10, 0)));
        panel.add(ruLabel);

        panel.setMaximumSize(panel.getPreferredSize());

        return panel;
    }


    JPanel createInputFilesPanel() {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.X_AXIS));
        br1 = BorderFactory.createTitledBorder("");
        addBorder(br1, panel);


        JPanel buttonsPanel = new JPanel();
        buttonsPanel.setLayout(new BoxLayout(buttonsPanel, BoxLayout.Y_AXIS));
        panel.add(buttonsPanel);

        addFileButton = new JButton();
        changeFileButton = new JButton();
        removeFileButton = new JButton();
        addFileButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        changeFileButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        removeFileButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        addFileButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                addFilesAction();
            }
        });
        changeFileButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                changeFileAction();
            }
        });
        removeFileButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                removeFilesAction();
            }
        });
        buttonsPanel.add(addFileButton);
        buttonsPanel.add(changeFileButton);
        buttonsPanel.add(removeFileButton);


        inFilesTable = new JTable(new AbstractTableModel() {
            public String getColumnName(int col) { return ""; }
            public int getRowCount() {
                if (inputFiles.get() == null) {
                    return 0;
                }
                return inputFiles.get().length;
            }
            public int getColumnCount() { return 1; }
            public Object getValueAt(int row, int col) {
                return inputFiles.get()[row];
            }
            public boolean isCellEditable(int row, int col) { return false; }
            public void setValueAt(Object value, int row, int col) {}
        });
        inFilesTable.setTableHeader(null);
        JScrollPane scrollPane = new JScrollPane(inFilesTable);
        scrollPane.setPreferredSize(new Dimension(800, 100));
        panel.add(scrollPane);


        return panel;
    }

    JPanel createParamsPanel() {
        JPanel panel = new JPanel();
        panel.setLayout(new FlowLayout());
        br2 = BorderFactory.createTitledBorder("");
        addBorder(br2, panel);


        final JPanel paramsPanel = new JPanel();
        paramsPanel.setLayout(new BoxLayout(paramsPanel, BoxLayout.Y_AXIS));
        paramsPanel.setBorder(BorderFactory.createEtchedBorder());
//        paramsPanel.setPreferredSize(new Dimension(360, 105));
        panel.add(paramsPanel);

        final JPanel optionsPanel = new JPanel();
        optionsPanel.setLayout(new BoxLayout(optionsPanel, BoxLayout.Y_AXIS));
        optionsPanel.setBorder(BorderFactory.createEtchedBorder());
//        optionsPanel.setPreferredSize(new Dimension(600, 115));
        panel.add(optionsPanel);


        // adding content to params panel
        kmerLabel = new JLabel();
        kmerSpinner = new JSpinner(new SpinnerNumberModel(31, 1, 31, 1));
        final JPanel kmerPanel = createSpinnerPanel(kmerLabel, kmerSpinner, 50);
        paramsPanel.add(kmerPanel);

        maxBadFrequencyLabel = new JLabel();
        maxBadFrequencySpinner = new JSpinner(new SpinnerNumberModel(1, 0, 100, 1));
        final JPanel maxBadFrequencyPanel = createSpinnerPanel(maxBadFrequencyLabel, maxBadFrequencySpinner, 50);

        minSeqLenLabel = new JLabel();
        minSeqLenSpinner = new JSpinner(new SpinnerNumberModel(100, 0, 10000, 1));;
        final JPanel minSeqLenPanel = createSpinnerPanel(minSeqLenLabel, minSeqLenSpinner, 100);

        minCompSizeLabel = new JLabel();
        minCompSizeSpinner = new JSpinner(new SpinnerNumberModel(1000, 0, 1000000, 1));
        final JPanel minCompSizePanel = createSpinnerPanel(minCompSizeLabel, minCompSizeSpinner, 100);

        maxCompSizeLabel = new JLabel();
        maxCompSizeSpinner = new JSpinner(new SpinnerNumberModel(10000, 0, 10000000, 1));
        final JPanel maxCompSizePanel = createSpinnerPanel(maxCompSizeLabel, maxCompSizeSpinner, 100);

        wrCheckBox = new JCheckBox();
        final JPanel wrPanel = createCheckBoxPanel(wrCheckBox);



        final JPanel more1Panel = new JPanel();
        more1Panel.setLayout(new FlowLayout(FlowLayout.CENTER));
        more1Button = new JButton();
        more1Button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (more1State == 0) {
                    paramsPanel.remove(more1Panel);

                    paramsPanel.add(maxBadFrequencyPanel);
                    paramsPanel.add(minSeqLenPanel);
                    paramsPanel.add(minCompSizePanel);
                    paramsPanel.add(maxCompSizePanel);
                    paramsPanel.add(wrPanel);

                    paramsPanel.add(more1Panel);
                    more1State = 1;
                    updateUIText();
                } else {
                    paramsPanel.remove(maxBadFrequencyPanel);
                    paramsPanel.remove(minSeqLenPanel);
                    paramsPanel.remove(minCompSizePanel);
                    paramsPanel.remove(maxCompSizePanel);
                    paramsPanel.remove(wrPanel);

                    more1State = 0;
                    updateUIText();
                }
            }
        });
        more1Panel.add(more1Button);
        paramsPanel.add(more1Panel);


        // adding content to options panel
        workDirLabel = new JLabel();
        JPanel workDirPanel = createLinePanel(workDirLabel);
        workDirField = new JTextField(20);
        workDirField.setEditable(true);
        setWorkDirButton = new JButton();
        setWorkDirButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                setWorkDirAction();
            }
        });
        workDirPanel.add(workDirField);
        workDirPanel.add(setWorkDirButton);
        optionsPanel.add(workDirPanel);

        memoryLabel = new JLabel();
        final JPanel memoryPanel = createLinePanel(memoryLabel);
        memoryField = new JTextField(7);
        setMemoryButton = new JButton();
        setMemoryButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                setMemoryAction();
            }
        });
        memoryPanel.add(memoryField);
        memoryPanel.add(setMemoryButton);

        avProcLabel = new JLabel();
        avProcSpinner = new JSpinner(new SpinnerNumberModel(1, 1, Runtime.getRuntime().availableProcessors(), 1));
        final JPanel avProcPanel = createSpinnerPanel(avProcLabel, avProcSpinner, 50);

        eaCheckBox = new JCheckBox();
        final JPanel eaPanel = createCheckBoxPanel(eaCheckBox);

        verboseCheckBox = new JCheckBox();
        final JPanel verbosePanel = createCheckBoxPanel(verboseCheckBox);
        verboseCheckBox.setPreferredSize(new Dimension(400, 20));

        final JPanel more2Panel = new JPanel();
        more2Panel.setLayout(new FlowLayout(FlowLayout.CENTER));
        more2Button = new JButton();
        more2Button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (more2State == 0) {
                    optionsPanel.remove(more2Panel);

                    optionsPanel.add(memoryPanel);
                    optionsPanel.add(avProcPanel);
                    optionsPanel.add(eaPanel);
                    optionsPanel.add(verbosePanel);

                    optionsPanel.add(more2Panel);
                    more2State = 1;
                    updateUIText();
                } else {
                    optionsPanel.remove(memoryPanel);
                    optionsPanel.remove(avProcPanel);
                    optionsPanel.remove(eaPanel);
                    optionsPanel.remove(verbosePanel);

                    more2State = 0;
                    updateUIText();
                }
            }
        });
        more2Panel.add(more2Button);
        optionsPanel.add(more2Panel);

        return panel;
    }

    JPanel createControlPanel() {
        JPanel panel = new JPanel();
        br3 = BorderFactory.createTitledBorder("");
        addBorder(br3, panel);


        runButton = new JButton();
        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                runAction();
            }
        });
        stopButton = new JButton();
        stopButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                stopAction();
            }
        });

        panel.add(runButton);
        panel.add(stopButton);
        
        panel.setMaximumSize(new Dimension(2000, (int) Math.ceil(panel.getPreferredSize().getHeight())));

        return panel;
    }

    JPanel createStatusPanel() {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        br4 = BorderFactory.createTitledBorder("");
        addBorder(br4, panel);


        runningStageText = new JLabel();
        runningStage = new JLabel();
        JPanel runningStatePanel = createStatusLine(runningStageText, runningStage);
        panel.add(runningStatePanel);

        elapsedTimeText = new JLabel();
        elapsedTime = new JLabel();
        JPanel elapsedTimePanel = createStatusLine(elapsedTimeText, elapsedTime);
        panel.add(elapsedTimePanel);

        remainingTimeText = new JLabel();
        remainingTime = new JLabel();
        JPanel remainingTimePanel = createStatusLine(remainingTimeText, remainingTime);
        panel.add(remainingTimePanel);

        overallProgressText = new JLabel();
        overallProgress = new JLabel();
        JPanel overallProgressPanel = createStatusLine(overallProgressText, overallProgress);
        panel.add(overallProgressPanel);
        
        logText = new JLabel();
        JPanel logPanel = createStatusLine(logText, null);
        logText.setPreferredSize(new Dimension(600, 20));
        log = new JTextArea(5, 40);
        log.setEditable(false);
        log.setBackground(panel.getBackground());
        JScrollPane logScrollPane = new JScrollPane(log);

        panel.add(logPanel);
        panel.add(logScrollPane);

        return panel;
    }



    // ============================ UI updating methods ===================================

    Dimension lastSize = null;
    Dimension newSize = null;

    void updateUIText() {
        // borders
        br1.setTitle(UIText[0][curLang]);
        br2.setTitle(UIText[1][curLang]);
        br3.setTitle(UIText[2][curLang]);
        br4.setTitle(UIText[3][curLang]);

        // buttons
        addFileButton.setText(UIText[4][curLang]);
        changeFileButton.setText(UIText[5][curLang]);
        removeFileButton.setText(UIText[6][curLang]);
        more1Button.setText(UIText[7 + more1State][curLang]);
        more2Button.setText(UIText[7 + more2State][curLang]);
        setWorkDirButton.setText(UIText[9][curLang]);
        setMemoryButton.setText(UIText[10][curLang]);
        runButton.setText(UIText[11][curLang]);
        stopButton.setText(UIText[12][curLang]);

        // labels and checkboxes
        setLabelText(kmerLabel, kmerSize);
        setLabelText(maxBadFrequencyLabel, maxBadFrequency);
        setLabelText(minSeqLenLabel, minSequenceLength);
        setLabelText(minCompSizeLabel, minComponentSize);
        setLabelText(maxCompSizeLabel, maxComponentSize);
        setCheckBoxText(wrCheckBox, withoutRenumbering);

        setLabelText(workDirLabel, workDir);
        memoryLabel.setText(UIOptionsText[0][curLang]);
        setLabelText(avProcLabel, availableProcessors);
        eaCheckBox.setText(UIOptionsText[1][curLang]);
        setCheckBoxText(verboseCheckBox, verbose);
        ToolTipManager.sharedInstance().setInitialDelay(500);
        ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE);


        if (frame != null) {
            Dimension curSize = frame.getSize();
            Dimension prefSize = frame.getPreferredSize();
            Dimension size = new Dimension((int) Math.ceil(Math.max(prefSize.getWidth(), curSize.getWidth())),
                    (int) Math.ceil(Math.max(prefSize.getHeight(), curSize.getHeight())));
            
            if (lastSize != null && !newSize.equals(curSize)) {
                // user change the size, no need to restore last size
                lastSize = null;
            }

            if (!size.equals(curSize)) {
                // increasing size
                if (lastSize == null) {
                    lastSize = curSize;
                }
                newSize = size;
            } else {
                if (lastSize == null) {
                    // nothing to do
                } else {
                    if (lastSize.getHeight() >= prefSize.getHeight() && lastSize.getWidth() >= prefSize.getWidth()) {
                        size = lastSize;
                        lastSize = null;
                    } else {
                        size = new Dimension((int) Math.ceil(Math.max(prefSize.getWidth(), lastSize.getWidth())),
                                (int) Math.ceil(Math.max(prefSize.getHeight(), lastSize.getHeight())));
                        newSize = size;
                    }
                }
            }

            frame.setSize(size);
            frame.validate();
//            frame.repaint();
        }
        
        // status
        runningStageText.setText(UIStatusText[0][curLang]);
        elapsedTimeText.setText(UIStatusText[1][curLang]);
        remainingTimeText.setText(UIStatusText[2][curLang]);
        overallProgressText.setText(UIStatusText[3][curLang]);
        logText.setText(UIStatusText[4][curLang] + " (" + metafast.workDir.get() + File.separator + "log" + "):");
        updateStatusInfo();
    }
    
    void setLabelText(JLabel label, Parameter p) {
        if (curLang == 0) {
            label.setText(p.description.descriptionShort);
            label.setToolTipText(p.description.description);
        } else {
            label.setText(p.description.descriptionRuShort);
            label.setToolTipText(p.description.descriptionRu);
        }
    }

    void setCheckBoxText(JCheckBox checkBox, Parameter p) {
        if (curLang == 0) {
            checkBox.setText(p.description.descriptionShort);
            checkBox.setToolTipText(p.description.description);
        } else {
            checkBox.setText(p.description.descriptionRuShort);
            checkBox.setToolTipText(p.description.descriptionRu);
        }
    }

    void updateUIValuesFromGlobalValues() {
        kmerSpinner.setValue(kmerSize.get());
        maxBadFrequencySpinner.setValue(maxBadFrequency.get());
        minSeqLenSpinner.setValue(minSequenceLength.get());
        minCompSizeSpinner.setValue(minComponentSize.get());
        maxCompSizeSpinner.setValue(maxComponentSize.get());
        wrCheckBox.setSelected(withoutRenumbering.get());

        workDirField.setText(workDir.get());
        memoryField.setText(memory);
        avProcSpinner.setValue(availableProcessors.get());
        eaCheckBox.setSelected(ea);
        verboseCheckBox.setSelected(verbose.get());
    }

    void updateGlobalValuesFromUIValues() {
        kmerSize.set((Integer) kmerSpinner.getValue());
        maxBadFrequency.set((Integer) maxBadFrequencySpinner.getValue());
        minSequenceLength.set((Integer) minSeqLenSpinner.getValue());
        minComponentSize.set((Integer) minCompSizeSpinner.getValue());
        maxComponentSize.set((Integer) maxCompSizeSpinner.getValue());
        withoutRenumbering.set(wrCheckBox.isSelected());

        workDir.set(workDirField.getText());
        memory = memoryField.getText();
        availableProcessors.set((Integer) avProcSpinner.getValue());
        ea = eaCheckBox.isSelected();
        verbose.set(verboseCheckBox.isSelected());
    }



    // =========================== creating panels methods ================================

    void addBorder(TitledBorder b, JPanel panel) {
        panel.setBorder(
                BorderFactory.createCompoundBorder(
                        BorderFactory.createEmptyBorder(10,5,0,5),
                        BorderFactory.createCompoundBorder(
                                b,
                                BorderFactory.createEmptyBorder(1,3,1,3)
                        )
                ));
    }

    JPanel createLinePanel(JLabel label) {
        JPanel panel = new JPanel();
        panel.setLayout(new FlowLayout(FlowLayout.LEFT));
        label.setPreferredSize(new Dimension(250, 20));
        panel.add(label);
        return panel;
    }

    JPanel createSpinnerPanel(JLabel label, JSpinner spinner, int preferredWidth) {
        spinner.setPreferredSize(new Dimension(preferredWidth, 20));
        return createSpinnerPanel(label,spinner);
    }
    JPanel createSpinnerPanel(JLabel label, JSpinner spinner) {
        JPanel panel = createLinePanel(label);
        panel.add(spinner);
        return panel;
    }

    JPanel createCheckBoxPanel(JCheckBox checkBox) {
        JPanel panel = new JPanel();
        panel.setLayout(new FlowLayout(FlowLayout.LEFT));
        checkBox.setPreferredSize(new Dimension(300, 20));
        panel.add(checkBox);
        return panel;
    }

    JPanel createStatusLine(JLabel textLabel, JLabel statusLabel) {
        JPanel panel = new JPanel();
        panel.setLayout(new FlowLayout(FlowLayout.LEFT));

        textLabel.setPreferredSize(new Dimension(250, 20));

        panel.add(Box.createRigidArea(new Dimension(10, 0)));
        panel.add(textLabel);
        if (statusLabel != null) {
            panel.add(Box.createRigidArea(new Dimension(20, 0)));
            panel.add(statusLabel);
        }
        
        panel.setMaximumSize(new Dimension(2000, 20));

        return panel;
    }


    // ========================== buttons action methods ========================

    void addFilesAction() {
        JFileChooser fc = new JFileChooser(new File("."));
        fc.setMultiSelectionEnabled(true);
        int res = fc.showOpenDialog(GUI.this);
        if (res == JFileChooser.APPROVE_OPTION) {
            File[] oldFiles = inputFiles.get();
            if (oldFiles == null) {
                oldFiles = new File[]{};
            }
            File[] newFiles = fc.getSelectedFiles();
            File[] allFiles = new File[oldFiles.length + newFiles.length];
            System.arraycopy(oldFiles, 0, allFiles, 0, oldFiles.length);
            System.arraycopy(newFiles, 0, allFiles, oldFiles.length, newFiles.length);

            inputFiles.set(allFiles);
            inFilesTable.updateUI();
            inFilesTable.clearSelection();
        }
    }

    void changeFileAction() {
        if (inFilesTable.getSelectedRowCount() != 1) {
            JOptionPane.showMessageDialog(frame, UIText[25][curLang]);
        } else {
            File f = inputFiles.get()[inFilesTable.getSelectedRow()];
            JFileChooser fc = new JFileChooser(f);
            fc.setMultiSelectionEnabled(false);
            fc.setSelectedFile(f);
            int res = fc.showOpenDialog(GUI.this);
            if (res == JFileChooser.APPROVE_OPTION) {
                inputFiles.get()[inFilesTable.getSelectedRow()] = fc.getSelectedFile();
                inFilesTable.updateUI();
                inFilesTable.clearSelection();
            }
        }
    }

    void removeFilesAction() {
        if (inFilesTable.getSelectedRowCount() == 0) {
            JOptionPane.showMessageDialog(frame, UIText[26][curLang]);
        } else {
            File[] files = inputFiles.get();
            File[] newFiles = new File[files.length - inFilesTable.getSelectedRowCount()];
            int j = 0;
            for (int i = 0; i < files.length; i++) {
                if (!inFilesTable.isRowSelected(i)) {
                    newFiles[j] = files[i];
                    j++;
                }
            }
            inputFiles.set(newFiles);
            inFilesTable.updateUI();
            inFilesTable.clearSelection();
        }
    }

    void setWorkDirAction() {
        JFileChooser fc = new JFileChooser(workDir.get());
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        int res = fc.showOpenDialog(GUI.this);
        if (res == JFileChooser.APPROVE_OPTION) {
            workDir.set(fc.getSelectedFile().toString());
            workDirField.setText(workDir.get());
        }
    }



    // ================================ main methods ===================================

    void setMemoryAction() {
        memory = memoryField.getText();

        // trying to start jvm with such memory
        int exitValue = -1;
        Process p = null;

        try {
            String command = "java -Xmx" + memory + " -Xms" + memory + " -version";
            p = Runtime.getRuntime().exec(command);
            exitValue = p.waitFor();
        } catch (InterruptedException e1) {
        } catch (IOException e1) {
        }
        if (exitValue != 0) {
            String erS = TextUtils.getLogText(p);
            JOptionPane.showMessageDialog(frame, UIText[14][curLang] + NL + erS,
                    "Invalid memory size", JOptionPane.ERROR_MESSAGE);
            return;
        }

        if (runningThread != null) {
            JOptionPane.showMessageDialog(frame, UIText[15][curLang],
                    "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        restartWithNewJVMOptions();
    }

    void restartWithNewJVMOptions() {
        assert runningThread == null;
        frame.setVisible(false);
        updateGlobalValuesFromUIValues();


        // generating command to run
        String metafastPath = Runner.getJarFilePath();

        jvmArgsStr += "-Xmx" + memory + " -Xms" + memory + " ";
        if (ea) {
            jvmArgsStr += "-ea ";
        }

        StringBuilder sb = new StringBuilder();
        for (Parameter p : almostAllParams) {
            if (p.description.tClass == Boolean.class) {
                if (p.get() != null && ((Boolean) p.get())) {
                    sb.append("--" + p.description.name + " ");
                }
            } else {
                String vs = Tool.objectToString(p.get());
                if (vs.contains(" ") && (vs.charAt(0) != '"' || vs.charAt(vs.length()-1) != '"')) {
                    vs = '"' + vs + '"';
                }
                sb.append("--" + p.description.name + " " + vs + " ");
            }
        }
        if (inputFiles.get() != null) {
            sb.append("--" + inputFiles.description.name + " " + TextUtils.arrayToString(inputFiles.get(), " "));
        }
        String args = sb.toString();


        String command = "java " + jvmArgsStr + "-jar " + metafastPath + " -gui " + args;
        System.err.println("command = " + command);
        
        try {
            Process p = Runtime.getRuntime().exec(command);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.exit(0);
    }


    // =================================== launching methods ==================================

    Thread runningThread, serviceThread;

    AppenderSkeleton guiApp = new AppenderSkeleton() {
        PatternLayout layout = new PatternLayout("%d{yyyy.MM.dd HH:mm:ss} %p %c: %m%n");
        @Override
        protected void append(LoggingEvent loggingEvent) {
            log.append(layout.format(loggingEvent));
        }
        @Override
        public void close() {
        }
        @Override
        public boolean requiresLayout() {
            return false;
        }
    };
    {
        Logger.getRootLogger().addAppender(guiApp);
    }


    
    void runAction() {
        if (runningThread != null) {
            JOptionPane.showMessageDialog(frame, UIText[16][curLang], "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        // preparing
        updateGlobalValuesFromUIValues();
        startNow.set(true);
        if (!memory.equals(currentMemory)) {
            setMemoryAction();
        }
        if (ea != currentEA) {
            restartWithNewJVMOptions();
        }
        startNow.set(false);

        for (Thread t : Thread.getAllStackTraces().keySet()) {
            threads.add(t.getName());
        }

        // updating global variables...
        Tool.launchedFromGUI = true;
        Tool.isInterrupted = false;
        Tool.lastException = null;

        Tool.forceParameter.set(false);
        Tool.continueParameter.set(false);
        File inputParamFile = new File(workDir.get(), Tool.IN_PARAM_FILE);
        if (inputParamFile.exists()) {
            Object[] options = {UIText[18][curLang], UIText[19][curLang], UIText[20][curLang]};
            int res = JOptionPane.showOptionDialog(frame, UIText[17][curLang].replace("WD", workDir.get() + File.separator), "Question",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.QUESTION_MESSAGE,
                    null,
                    options,
                    options[2]);
            if (res == 0) {
                Tool.forceParameter.set(true);
            } else if (res == 1) {
                Tool.continueParameter.set(true);
            } else {
                return;
            }
        }
        if (Tool.forceParameter.get() && (inputFiles.get() == null || inputFiles.get().length == 0)) {
            JOptionPane.showMessageDialog(frame, UIText[22][curLang], "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        if (verbose.get()) {
            guiApp.setThreshold(Level.DEBUG);
        } else {
            guiApp.setThreshold(Level.INFO);
        }
        
        logText.setText(UIStatusText[4][curLang] + " (" + metafast.workDir.get() + File.separator + "log" + "):");
        log.setText("");

        
        // launching ...
        metafast.progress.reset();
        updateStatusInfo();
        runningThread = new Thread(new Runnable() {
            @Override
            public void run() {
                metafast.mainImpl(new String[]{});
            }
        });
        runningThread.start();

        String serviceThreadName = "Service-Thread-1";
        serviceThread = new Thread(new Runnable() {
            @Override
            public void run() {
                while (runningThread.isAlive()) {
                    updateStatusInfo();
                    try {
                        Thread.sleep(1000);
                    } catch (InterruptedException e) {
                        System.err.println("Service thread was interupted.");
                    }
                }
                updateStatusInfo();
                runningThread = null;
                if (!Tool.isInterrupted && Tool.lastException == null) {    // finished normally
                    runningStage.setText("Finished!");
                    runningStage.setForeground(new Color(0,128,0));     // dark green
                    JOptionPane.showMessageDialog(frame,
                            UIText[24][curLang],
                            "Finished!", JOptionPane.INFORMATION_MESSAGE);
                } else {
                    updateStatusInfo();
                    if (Tool.lastException != null) {
                        JOptionPane.showMessageDialog(frame, UIText[23][curLang] + NL + Tool.lastException, "Error", JOptionPane.ERROR_MESSAGE);
                    }
                }
            }
        }, serviceThreadName);
        threads.add(serviceThreadName);
        serviceThread.start();
    }

    HashSet<String> threads = new HashSet<String>();

    void stopAction() {
        if (runningThread == null) {
            JOptionPane.showMessageDialog(frame, UIText[21][curLang], "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        Tool.isInterrupted = true;
        for (Thread t : Thread.getAllStackTraces().keySet()) {
            if (!threads.contains(t.getName())) {
                t.interrupt();
            }
        }
    }

    void updateStatusInfo() {
        if (metafast.progress.currentRunningStep != null) {
            String s = metafast.progress.currentRunningStep.getStageName(curLang) + " " +
                    "(" + ((curLang == 0) ? "stage " : "этап ") + (metafast.progress.doneSteps + 1) + ((curLang == 0) ? " of " : " из ") + metafast.progress.allSteps +")";
            if (Tool.isInterrupted) {
                if (runningThread != null) {
                    s += (curLang == 0) ? " (interrupting...)" : " (прерывание...)";
                } else {
                    s += (curLang == 0) ? " (interrupted)" : " (прерван)";
                    runningStage.setForeground(Color.RED);
                }
            }
            runningStage.setText(s);
        }

        if (runningThread == null && !Tool.isInterrupted) {
            // clearing the status
            runningStage.setText("");
            runningStage.setForeground(Color.BLACK);
            remainingTime.setText("");
            return;
        }

        elapsedTime.setText(metafast.progress.getRunningTime());
        remainingTime.setText(metafast.progress.getRemainingTime());
        overallProgress.setText(String.format("%.1f%%", metafast.progress.progress * 100.0));
    }


    
    
    // ========================== helpful methods =============================

    String jvmArgsStr = "";

    void determineMemoryDefaultValue() {
        currentMemory = null;
        RuntimeMXBean runtimeMXBean = ManagementFactory.getRuntimeMXBean();
        List<String> arguments = runtimeMXBean.getInputArguments();

        Pattern pattern1 = Pattern.compile("(?i)^-xmx(.*)$");
        Pattern pattern2 = Pattern.compile("(?i)^-xms(.*)$");
        Pattern pattern3 = Pattern.compile("(?i)^-ea$");
        for (String s : arguments) {
            Matcher matcher = pattern1.matcher(s);
            if (matcher.matches()) {
                currentMemory = matcher.group(1);
            } else {
                if (!pattern2.matcher(s).matches() && !pattern3.matcher(s).matches()) {
                    jvmArgsStr += s + " ";
                }
            }
        }
        if (currentMemory == null) {
            currentMemory = Misc.availableMemoryAsStringForJVM();
        }
        memory = currentMemory;
    }

    void determineEADefaultValue() {
        currentEA = false;
        assert testEA();
        ea = currentEA;
    }

    boolean testEA() {
        currentEA = true;
        return true;
    }


    static JFrame frame;

    public static void mainImpl(final String[] args) {
        Runner.getJarFilePath();  // checking that running from a jar file

        GUI gui = new GUI(args);

        frame = new JFrame("MetaFast (version " + Runner.getVersionNumber() + ")");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        frame.add(gui);

        frame.pack();
        frame.setVisible(true);

        if (gui.startNow.get()) {
            gui.runAction();
        }
    }
}
