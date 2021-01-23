/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vic;

import java.io.*;
import java.util.*;
import static java.util.Arrays.asList;
import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException;
import org.apache.commons.cli.*;

/**
 * VIC.java
 *
 * @author Yue Hu and Muqing Yan
 */
public class VIC {

    protected static Object key, value;
    protected static String inputfile, buildver, outputfile, inputtype, otherinfo;
    protected static String vicdb, database_humandb, evidencefile, cancer_type, knownlist;
    protected static int sum_sample, line_sum;
    protected static String lof_genes, tableannovar, convert_annovar, annotatevariation, mim2gene, mim_pheno, mim_orpha;
    protected static String orpha, exclude_snps, cgi_markers, add_markers;
    protected static String knowgenecanonical, cgiline2;
    protected static String cancer_pathways, cancer_genes, civic_markers, cancer_types, outfilename;

    protected static Map paras = new HashMap();

    // some important datasets/lists
    protected static Map user_evidence_dict = new HashMap();
    protected static Map lof_genes_dict = new HashMap();
    protected static Map mim2gene_dict = new HashMap();
    protected static Map mim2gene_dict2 = new HashMap();
    protected static Map exclude_snps_dict = new HashMap();
    protected static Map mim_pheno_dict = new HashMap();
    protected static Map mim_orpha_dict = new HashMap();
    protected static Map orpha_dict = new HashMap();
    protected static Map knownGeneCanonical_dict = new HashMap();
    protected static Map knownGeneCanonical_st_dict = new HashMap();
    protected static Map knownGeneCanonical_ed_dict = new HashMap();
    protected static Map add_markers_dict = new HashMap();
    protected static Map cgi_markers_dict = new HashMap();
    protected static Map civ_markers_dict = new HashMap();
    protected static Map cancer_pathways_dict = new HashMap();
    protected static Map cancer_genes_dict = new HashMap();
    protected static Map cancer_types_dict = new HashMap();

    protected static Map<String, Integer> Freqs_flgs = new HashMap<>();
    protected static Map<String, Integer> Funcanno_flgs = new HashMap<>();
    protected static Map<String, Integer> Allels_flgs = new HashMap<>();

    protected static String[] add_d, civ_d;
    protected static String[] cgi_d;
    protected static List<String> known_list = new ArrayList<>();
    protected static String cosmicVersion = "cosmic84_coding";
    protected static String dbnsfpVersion = "dbnsfp35a";
    protected static String clinvarVersion = "clinvar_20190305";
    protected static String gnomadVersion = "gnomad211_exome";

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     * @throws java.lang.InterruptedException
     * @throws org.apache.commons.cli.ParseException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, ParseException,
            InterruptedException {

        System.out.println("Notice: Your command of VIC is " + Arrays.toString(args));

        Options options = new Options();

        // VIC 
        options.addOption("h", "help", false, "help message");
        options.addOption("b", "buildver", true, "The genomic build version, it can be hg19 and will support GRCh37 hg18 GRCh38 later");
        options.addOption("i", "inputfile", true, "--required argument: The input file contains your variant");
        options.addOption("o", "outputfile", true, "--required argument: The prefix of output file which contains the results, the file of results will be as [$$prefix].vic");
        options.addOption("input_type", true, "--required argument: The input file type, it can be AVinput(Annovar's format), VCF(VCF with single sample), vcf_m(VCF with multiple samples)");

        // VIC Other Options
        options.addOption("db", "database_vic", true, "--required argument: <path>The database location/dir for the VIC dataset file");
        options.addOption("s", "evidence_file", true, "User specified Evidence file for each variant");
        options.addOption("l", "user_specified_list", true, "User specified variant list with self-designated clinical impact significance, must be in AVinput format");
        options.addOption("cancer_type", true, "Please check the 'vic.cancer.types' in vicdb for available cancer types. All the cancer types should be in the shout-cut form, e.g.:AA. Or include cancer types in the otherinfo column of the AVinput");

        options.addOption("cosmic", "cosmic_version", true, "User specified cosmic version");
        options.addOption("dbnsfp", "dbnsfp_version", true, "User specified dbnsfp version");
        options.addOption("clinvar", "clinvar_version", true, "User specified clinvar version");
        options.addOption("gnomad", "gnomad_version", true, "User specified gnomad version");

        // Annovar Options, Caution: check these options from manual of Annovar
        options.addOption("table_annovar", true, "<path>The Annovar perl script of table_annovar.pl");
        options.addOption("convert2annovar", true, "<path>The Annovar perl script of convert2annovar.pl");
        options.addOption("annotate_variation", true, "<path>The Annovar perl script of annotate_variation.pl");
        options.addOption("d", "database_locat/humandb", true, "<path>The database location/dir for the annotation dataset");

        options.addOption("skip_annovar", false, "Skip the Annovar annotation, this can be true only after you already got the annovar's annotation result");
        options.addOption("otherinfo", true, "true or false:print out otherinfo (infomration in fifth column in queryfile, default: TRUE). "
                + "this option only perform well with AVinput file, and the other information only can be put in the fifth column. The information in >5th column will be lost."
                + " When input as VCF or VCF_m files with otherinfo option, only het/hom will be kept, depth and qual will be lost, the cancer type should be provide by command option.");

// System.out.println(Arrays.toString(args));
        if (Arrays.toString(args).equals("[]")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("java -jar VIC.jar", options, true);
            System.exit(0);
        }

// CommandLineParser parser = new PosixParser();
        CommandLineParser parser = new DefaultParser();
        CommandLine option = parser.parse(options, args);

// System.out.println(options.toString());
// if (options.hasOption("")) {
// HelpFormatter formatter = new HelpFormatter();
// formatter.printHelp("java -jar VIC.jar", options, true);
// System.exit(0);
// }
//
        if (option.hasOption("h")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("java -jar VIC.jar", options, true);
            System.exit(0);
        }

        if (option.hasOption("i")) {
            inputfile = option.getOptionValue("i");
            paras.put("inputfile", option.getOptionValue("i"));
        }

        if (option.hasOption("o")) {
            outputfile = option.getOptionValue("o");
            paras.put("outputfile", option.getOptionValue("o"));
        }

        if (option.hasOption("input_type")) {
            inputtype = option.getOptionValue("input_type");
            paras.put("input_type", option.getOptionValue("input_type"));
        }

        if (!option.hasOption("cancer_type")) {
            cancer_type = "CANCER";
        } else {
            cancer_type = option.getOptionValue("cancer_type");
            paras.put("cancer_type", cancer_type);
        }

        if (option.hasOption("cosmic")) {
            cosmicVersion = option.getOptionValue("cosmic_version");
            if (cosmicVersion.lastIndexOf("_coding") < 0) {
                cosmicVersion += "_coding";
            }
        } else {
            cosmicVersion = "cosmic84_coding";
        }

        if (option.hasOption("dbnsfp")) {
            dbnsfpVersion = option.getOptionValue("dbnsfp_version");
        } else {
            dbnsfpVersion = "dbnsfp35a";
        }

        if (option.hasOption("clinvar")) {
            clinvarVersion = option.getOptionValue("clinvar_version");
        } else {
            clinvarVersion = "clinvar_20190305";
        }

        if (option.hasOption("gnomad")) {
            gnomadVersion = option.getOptionValue("gnomad_version");
        } else {
            gnomadVersion = "gnomad211_exome";
        }

        if (option.hasOption("b")) {
            buildver = option.getOptionValue("b");
        } else {
            buildver = "hg19";
        }
        paras.put("buildver", buildver);

        if (option.hasOption("d")) {
            database_humandb = option.getOptionValue("d");
        } else {
            database_humandb = "./humandb";
        }
        paras.put("humandb", database_humandb);

        if (option.hasOption("db")) {
            vicdb = option.getOptionValue("db");
            lof_genes = vicdb + "/LOF.genes.exac_me_cancers";
            mim2gene = vicdb + "/mim2gene.txt";
            mim_pheno = vicdb + "/mim_pheno.txt";
            mim_orpha = vicdb + "/mim_orpha.txt";
            orpha = vicdb + "/orpha.txt";
            knowgenecanonical = vicdb + "/knownGeneCanonical.txt";
            exclude_snps = vicdb + "/ext.variants";
            cgi_markers = vicdb + "/cgi_biomarkers_20180117.txt";
            add_markers = vicdb + "/add_marker.full.txt";
            cancer_pathways = vicdb + "/cancer_pathways.list";
            cancer_genes = vicdb + "/cancer_genes.list";
            civic_markers = vicdb + "/civic_2019_06.txt";//BY MQY
            cancer_types = vicdb + "/vic.cancer.types";
        } else {
            vicdb = "./vicdb";
            lof_genes = vicdb + "/LOF.genes.exac_me_cancers";
            mim2gene = vicdb + "/mim2gene.txt";
            mim_pheno = vicdb + "/mim_pheno.txt";
            mim_orpha = vicdb + "/mim_orpha.txt";
            orpha = vicdb + "/orpha.txt";
            knowgenecanonical = vicdb + "/knownGeneCanonical.txt";
            exclude_snps = vicdb + "/ext.variants";
            cgi_markers = vicdb + "/cgi_biomarkers_20180117.txt";
            add_markers = vicdb + "/add_marker.full.txt";
            cancer_pathways = vicdb + "/cancer_pathways.list";
            cancer_genes = vicdb + "/cancer_genes.list";
            civic_markers = vicdb + "/civic_2019_06.txt";
            cancer_types = vicdb + "/vic.cancer.types";
        }
        paras.put("vicdb", vicdb);

        exclude_snps = exclude_snps + "." + buildver;

        // 09072018
        if (option.hasOption("s")) {
            evidencefile = option.getOptionValue("s");
            File evidence_file = new File(evidencefile);
            if (!evidence_file.isFile()) {
                System.out.println("Warning: Your specified evidence file " + evidence_file
                        + " is not here, please check the path of your evidence file.");
                System.out.println("Your analysis will begin without your specified evidence.");
            } else {
                System.out.println("INFO: Your analysis will begin with your provided evidence");
                String line2;
                try (BufferedReader reader_evi = new BufferedReader(new FileReader(evidence_file))) {
                    while ((line2 = reader_evi.readLine()) != null) {
                        String[] cls2 = line2.split("\t");
                        if (cls2.length > 1) {
                            String keys = cls2[0] + "_" + cls2[1] + "_" + cls2[2] + "_" + cls2[3];
                            keys = keys.replaceAll("(?i)chr", "");
// System.out.println("====the evidence key is " + keys);
                            user_evidence_dict.put(keys, cls2[4].toUpperCase());
                        }
                    }
                } catch (IOException e) {
                    System.out.println(e);
                }
            }
        } else {
            evidencefile = " ";
        }

        // 03142019 by Muqing Yan
        if (option.hasOption("l")) {
            knownlist = option.getOptionValue("l");
            File known_file = new File(knownlist);
            if (!known_file.isFile()) {
                System.out.println("Warning: Your specified known list file " + known_file
                        + " is not here, please check the path of your known list file.");
                System.out.println("Your analysis will begin without your specified evidence.");
                knownlist = " ";
            } else {
                System.out.println("==Waring: Your analysis will begin with your specified significance.");
                String line3;
                try (BufferedReader reader_list = new BufferedReader(new FileReader(known_file))) {
                    while ((line3 = reader_list.readLine()) != null) {
                        String[] cls2 = line3.split("\t");
                        if (cls2.length > 1) {
                            String known_variant = cls2[0] + "#" + cls2[1] + "#" + cls2[2] + "#" + cls2[3] + "#" + cls2[4] + "#" + cls2[5];
                            known_variant = known_variant.replaceAll("(?i)chr", "");
                            known_list.add(known_variant);
                        }
                    }
                } catch (IOException e) {
                    System.out.println(e);
                }
            }
        } else {
            knownlist = " ";
        }

        if (option.hasOption("table_annovar")) {
            File f1 = new File(option.getOptionValue("table_annovar"));
            if (f1.isFile()) {
                tableannovar = option.getOptionValue("table_annovar");

            } else {
                System.out.println("Error:The Annovar file " + tableannovar + " is not here, please download ANNOVAR firstly:http://www.openbioinformatics.org/annovar");
                System.exit(0);
            }
        }
        paras.put("table_annovar", option.getOptionValue("table_annovar"));

        if (!option.hasOption("otherinfo")) {
            paras.put("otherinfo", "TRUE");
        } else {
            paras.put("otherinfo", option.getOptionValue("otherinfo"));
        }

        if (option.hasOption("convert2annovar")) {
            File f2 = new File(option.getOptionValue("convert2annovar"));
            if (f2.isFile()) {
                convert_annovar = option.getOptionValue("convert2annovar");
            } else {
                convert_annovar = option.getOptionValue("convert2annovar");
                System.out.println("Error: The Annovar file" + convert_annovar + " is not here, please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar");
                System.exit(0);
            }
        }
        paras.put("convert2annovar", option.getOptionValue("convert2annovar"));

        if (option.hasOption("annotate_variation")) {
            File f3 = new File(option.getOptionValue("annotate_variation"));
            if (f3.isFile()) {
                annotatevariation = option.getOptionValue("annotate_variation");
            } else {
                System.out.println("Error: The Annovar file " + annotatevariation + " is not here, please download ANNOVAR firstlliney: http://www.openbioinformatics.org/annovar");
                System.exit(0);
            }
            paras.put("annotate_variation", option.getOptionValue("annotate_variation"));
        }

        File input_file = new File(inputfile);
        File evidence_file = new File(evidencefile);

        if (!input_file.isFile()) {
            System.out.println("Error: Your input file " + inputfile + " is not here, please check the path of your input file. ");
            System.exit(0);
        }
        if (option.hasOption("s") && !evidence_file.isFile()) {
            System.out.println("INFO: Your specified evidence file " + evidencefile + " is not here, please check the path of your evidence file.");
            System.out.println("Your analysis will begin without your specified evidence.");
        }

        System.out.println("INFO:The options are " + paras);

        Check_Downdb();
        Check_Input();

        if (!option.hasOption("skip_annovar")) {
            Check_annovar_result();
        } else {
            System.out.println("Warning: You activated the option of --skip_annovar, the Annovar will not run!");
            System.out.println("Warning: The VIC will interpret the variants based on your old annotation information!");
        }

        Read_dataset();

        int out_annf = 0;
        int o1 = -1;
        if (outputfile.contains("/")) {
            o1 = outputfile.lastIndexOf("/");
        }
        String oname = outputfile.substring(o1 + 1);
        String opath = outputfile.substring(0, o1 + 1);

        if (opath.equals("")) {
            opath = "./";
        }
// System.out.println("===#####oname and opath are: " + oname + "###" + opath);

        File output_path = new File(opath);

        String[] annovar_outfiles = output_path.list((File dir, String name) -> name.endsWith(buildver + "_multianno.txt") && name.startsWith(oname));

        for (String annovar_outfile : annovar_outfiles) {
            System.out.printf("annovar_outfile is %s%s \n", opath, annovar_outfile);
            int sum1 = Check_genes(opath + annovar_outfile);
            int sum2 = My_Cancer_var_can(opath + annovar_outfile);//process each ANNO file with VIC

            out_annf = out_annf + 1;
            outfilename = opath + annovar_outfile + ".vic";
            File outfile = new File(outfilename);
            System.out.println(outfile);
            if (outfile.isFile()) {
                System.out.printf("Notice: About %d lines in your variant file! \n", (sum1 - 1));
                System.out.printf("Notice: About %d variants has been processed by VIC \n", (sum2 - 1));
                if (!"vcf_m".equals(inputtype.toLowerCase())) {
                    System.out.println("Notice: The VIC is finished, the output file is [ " + outfilename + " ]");
                }
            }
        }

        if ("vcf_m".equals(inputtype.toLowerCase())) {
            System.out.println("Notice: The VIC for VCF with multiple samples is finished, the output files are as [ " + "annovar_outfile" + "<sample name>.vic ]:");
            sum_sample = 1;
            String[] f = output_path.list((File dir, String name) -> name.endsWith(buildver + "_multianno.txt.vic") && name.startsWith(oname));
            for (String f1 : f) {
                System.out.println(f1);
            }
        }

        if (out_annf == 0) {
            System.out.println("Warning: The VIC seems not run correctly, please check your inputs , options and configure file!");
            System.out.println("Error: The VIC did not find the annotation result file from ANNOVAR!");
            System.out.println("Error: The name of annotation result file should be like " + outputfile + ".<*sample_name*>." + buildver + "_multianno.txt");
        }

        System.out.println("Done.");
    }

    public static void Read_dataset() throws FileNotFoundException, IOException {
        // 0. read the user specified evidence file
        try {
            File evidence_file = new File(evidencefile);
            if (evidence_file.isFile()) {
                String line2;
                try (BufferedReader reader = new BufferedReader(new FileReader(evidence_file))) {
                    while ((line2 = reader.readLine()) != null) {
                        String[] cls2 = line2.split("\t");
                        if (cls2.length > 4) {
                            String keys = cls2[0] + "_" + cls2[1] + "_" + cls2[2] + "_" + cls2[3];
                            keys = keys.replaceAll("(?i)chr", "");
                            user_evidence_dict.put(keys, cls2[4].toUpperCase());
                        }
                    }
                } catch (IOException e) {

                }
            }
        } catch (Exception e) {
            System.out.println("Error:can't read the user sepcified evidence file" + evidencefile);
        }

        // 1.LOF gene list 
        try {

            File lof_gene = new File(lof_genes);
            String line2;
            try (BufferedReader reader_lof_gene = new BufferedReader(new FileReader(lof_gene))) {
                while ((line2 = reader_lof_gene.readLine()) != null) {
                    String[] cls2 = line2.split("\t");
                    // System.out.println(cls2[0]);
                    if (cls2[0].length() > 1) {
                        lof_genes_dict.put(cls2[0], "1");
                    }
                }
                reader_lof_gene.close();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Error: can\'t read the LOF genes file " + lof_genes);
            System.out.println("Error: Please download it from the source website");
            return;
        }

        // 2.cgi_markers
        try {
            String line2;
            File cgi_marker = new File(cgi_markers);
            try (BufferedReader reader_cgi_marker = new BufferedReader(new FileReader(cgi_marker))) {
                String default_s = "";
                int linecount = -1;
                String mutt, gene1, mutb;
                while ((line2 = reader_cgi_marker.readLine()) != null) {
                    cgi_d = line2.split("\t");
                    linecount += 1;
                    if (cgi_d[0].length() > 1) {
                        String gene = cgi_d[1];
                        String mut_type = cgi_d[2];
                        String mut = cgi_d[3];
                        String[] mut_types = mut_type.split(";");
                        int list_mut_size = mut_types.length;
                        String[] mut_alts = mut.split(";");
                        for (int i = 0; i < list_mut_size; i++) {
                            String mut_list = mut_alts[i];
                            String[] mut_lists = mut_list.split(":");
                            if (mut_lists.length > 1) {
                                gene1 = mut_lists[0];
                                mutb = mut_lists[1];
                            } else {
                                gene1 = gene;
                                mutb = mut_alts[i];
                            }
                            switch (mut_types[i]) {
                                case "FUS":
                                    mutt = "fus_" + mut_alts[i].replaceAll("__", "-");
                                    key = mutt;
                                    value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                    cgi_markers_dict.put(key, value);
                                    break;

                                case "EXPR":
                                    mutt = gene1 + "_expr_" + mutb;
                                    key = mutt;
                                    value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                    cgi_markers_dict.put(key, value);
                                    break;
                                case "CNA":
                                    mutt = gene1 + "_cna_" + mutb;
                                    key = mutt;
                                    value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                    cgi_markers_dict.put(key, value);
                                    break;
                                case "BIA":
                                    mutt = gene1 + "_bia_" + mutb;
                                    key = mutt;
                                    value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                    cgi_markers_dict.put(key, value);
                                    break;
                                case "MUT":
                                    for (String mutt_s : mutb.split(", ")) {
                                        mutt = gene1 + "_" + mutt_s;
                                        key = mutt;
                                        value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                        cgi_markers_dict.put(key, value);
                                    }
                                    break;
                                default:
                                    key = gene1 + "_" + mutb;
                                    value = linecount + ", " + cgi_markers_dict.put(key, default_s);
                                    cgi_markers_dict.put(key, value);
                                    break;
                            }
                        }
                    }
                }
                reader_cgi_marker.close();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Error: can\'t read the cgi_markers genes file " + cgi_markers);
            System.out.println("Error: Please download it from the source website");
        }

        // 3.add_markers from pmkb
        try {
            File add_marker = new File(add_markers);
            String line2;
            try (BufferedReader reader_add_marker = new BufferedReader(new FileReader(add_marker))) {
                int linecount = -1;
                String mut;
                String default_s;
                while ((line2 = reader_add_marker.readLine()) != null) {
                    linecount += 1;
                    add_d = line2.split("\t");
                    String gene = add_d[0];
                    mut = add_d[1].trim();
                    if (mut.toLowerCase().contains("any ")) {
                        if (mut.equalsIgnoreCase("any mutation")) {
                            mut = "mut_any";
                        } else if (mut.equalsIgnoreCase("any insertion")) {
                            mut = "insertion_any";
                        } else if (mut.equalsIgnoreCase("any deletion")) {
                            mut = "deletion_any";
                        } else if (mut.equalsIgnoreCase("any frameshift")) {
                            mut = "frameshift_any";
                        } else if (mut.equalsIgnoreCase("any indel")) {
                            mut = "indel_any";
                        } else if (mut.equalsIgnoreCase("any missense")) {
                            mut = "missenese_any";
                        } else if (mut.equalsIgnoreCase("any nonsense")) {
                            mut = "nonsense_any";
                        }
                        key = gene + "_" + mut;
                        default_s = "";
                        value = linecount + ", " + add_markers_dict.put(key, default_s);
                        add_markers_dict.put(key, value);
                    } else if (mut.toLowerCase().contains("rearrangement")) {
                        mut = "fus_" + gene;
                        key = mut;
                        default_s = "";
                        value = linecount + ", " + add_markers_dict.put(key, default_s);
                        add_markers_dict.put(key, value);
                    } else if (mut.toLowerCase().contains("copy number ")) {
                        String[] mutt = mut.split("copy number ");
                        String muts = "cna_" + mutt[1];
                        key = gene + "_" + muts;
                        default_s = "";
                        value = linecount + ", " + add_markers_dict.put(key, default_s);
                        add_markers_dict.put(key, value);
                    } else if (mut.toLowerCase().contains("exon(s)")) {
                        String[] str = mut.split("exon(s) | ");
                        if (!mut.toLowerCase().contains(", ")) {
                            String poss = str[1];
                            String mutb = str[2];
                            String muts = "exon_" + poss + "_" + mutb;
                            key = gene + "_" + muts;
                            default_s = " ";
                            value = linecount + ", " + add_markers_dict.put(key, default_s);
                            add_markers_dict.put(key, value);
                        } else {
                            str = mut.replace(", ", " ").replaceAll(" ", " ").split(" ");
                            int index1 = mut.lastIndexOf(", ");
                            int index2 = mut.lastIndexOf(" ");
                            String pos_end = mut.substring(index1 + 2, index2);
                            String mutb = mut.substring(index2 + 1);
                            String muts = "exon_" + pos_end + "_" + mutb;
                            key = gene + "_" + muts;
                            default_s = "";
                            value = linecount + ", " + add_markers_dict.put(key, default_s);
                            add_markers_dict.put(key, value);
// System.out.println(value);
                            for (int i = 1; i < str.length - 2; i++) {
                                String pos = str[i];
                                String muta = "exon_" + pos + "_" + mutb;
                                key = gene + "_" + muta;
                                default_s = " ";
                                value = linecount + ", " + add_markers_dict.put(key, default_s);
                                add_markers_dict.put(key, value);
                            }
                        }
                    } else if (mut.toLowerCase().contains("codon(s)")) {
                        String[] str = mut.split("codon(s) | ");
                        if (!mut.toLowerCase().contains(", ")) {
                            String poss = str[1];
                            String mutb = str[2];
                            String muts = "codon_" + poss + "_" + mutb;
                            key = gene + "_" + muts;
                            default_s = " ";
                            value = linecount + ", " + add_markers_dict.put(key, default_s);
                            add_markers_dict.put(key, value);
                        } else {
                            str = mut.replace(", ", " ").replaceAll(" ", " ").split(" ");
                            int index1 = mut.lastIndexOf(", ");
                            int index2 = mut.lastIndexOf(" ");
                            String pos_end = mut.substring(index1 + 2, index2);
                            String mutb = mut.substring(index2 + 1);
                            String muts = "codon_" + pos_end + "_" + mutb;
                            key = gene + "_" + muts;
                            default_s = " ";
                            value = linecount + ", " + add_markers_dict.put(key, default_s);
                            add_markers_dict.put(key, value);

                            for (int i = 1; i < str.length - 2; i++) {
                                String pos = str[i];
                                String muta = "codon_" + pos + "_" + mutb;
                                key = gene + "_" + muta;
                                default_s = " ";
                                value = linecount + ", " + add_markers_dict.put(key, default_s);
                                add_markers_dict.put(key, value);
                            }
                        }
                    } else {
                        key = gene + "_" + mut;
                        default_s = " ";
                        value = linecount + ", " + add_markers_dict.put(key, default_s);
                        add_markers_dict.put(key, value);
                    }
                }

                reader_add_marker.close();
            }
        } catch (FileNotFoundException e) {
            System.out.printf("Error: can\'t read the additional markers file %s", add_markers);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 4.civic_markers from civic// BY MQY
        try {
            File civic_marker = new File(civic_markers);
            try (BufferedReader reader_civic_marker = new BufferedReader(new FileReader(civic_marker))) {
                String line2;
                int linecount = -1;
                while ((line2 = reader_civic_marker.readLine()) != null) {
                    civ_d = line2.split("\t");
                    String gene = civ_d[0];
// System.out.println("=======civic gene: " + gene);
                    if (civ_d[0].length() > 1 && linecount >= 0) {
                        key = civ_d[23] + "#" + civ_d[24] + "#" + civ_d[25] + "#" + civ_d[26] + "#" + civ_d[27];
// System.out.println("=======civic key: " + key);
                        value = civ_d[0] + ", " + civ_d[2] + ", " + civ_d[3] + ", " + civ_d[6] + ", " + civ_d[8] + ", " + civ_d[9] + ", " + civ_d[10] + ", " + civ_d[11] + ", " + civ_d[12] + ", " + civ_d[16];
                        civ_markers_dict.put(key, value);

                    }
                    linecount += 1;
                }
                reader_civic_marker.close();
            }
        } catch (FileNotFoundException e) {

        }

        // 5.OMIM mim2gene.txt file
        try {
            File mim_2gene = new File(mim2gene);

            try (BufferedReader reader_mim2gene = new BufferedReader(new FileReader(mim_2gene))) {
                String line2;

                while ((line2 = reader_mim2gene.readLine()) != null) {
                    String[] cls2 = line2.split("\t");
                    String cls0;

                    if (cls2.length > 4) {
                        cls0 = cls2[4];
                    } else {
                        cls0 = "";
                    }

                    if (cls2.length > 1) {
                        String[] keys = cls0.split(", ");
                        key = keys[0];
                        value = cls2[0];
                        mim2gene_dict.put(key, value);
                    }

                    if (cls2.length > 3) {
                        key = cls2[3].toUpperCase();
                        value = cls2[0];
                        mim2gene_dict2.put(key, value);
                    }

                }

                reader_mim2gene.close();
            }

        } catch (FileNotFoundException e) {

        }

// // 6.read the user specified SNP list, the variants will pass the frequency check//???NO "ext.variants" file in the vicdb
// File excludesnps = new File(exclude_snps);
// if (excludesnps.isFile()) {
// try {
// try (BufferedReader reader_excludesnp = new BufferedReader(new FileReader(excludesnps))) {
// String line2;
// while ((line2 = reader_excludesnp.readLine()) != null) {
// String[] cls2 = line2.split("\t");
// if (cls2.length > 1) {
// String keys = cls2[0] + "_" + cls2[1] + "_" + cls2[2] + "_" + cls2[3];
// key = keys.replaceAll("(?i)chr", "");
// value = "1";
// exclude_snps_dict.put(key, value);
// }
// }
// }
//
// } catch (IOException e) {
// System.out.printf("Error: can\\'t read the user specified SNP list file %s", exclude_snps);
// }
// }
//
        // 7. knownGeneCanonical exon file # caution the buildver, now it is hg19
        try {
            File knowngene = new File(knowgenecanonical);
            try (BufferedReader reader_knowngene = new BufferedReader(new FileReader(knowngene))) {
                Map knowngene_dict = new HashMap();
                Map knowngene_st_dict = new HashMap();
                Map knowngene_ed_dict = new HashMap();
                String line2;
                while ((line2 = reader_knowngene.readLine()) != null) {
                    String[] cls2 = line2.split(" ");
                    if (cls2.length > 1) {
                        key = cls2[0];
                        knowngene_dict.put(key, cls2[1]);
                        knowngene_st_dict.put(key, cls2[2]);
                        knowngene_ed_dict.put(key, cls2[3]);
                    }
                }
                reader_knowngene.close();
            }

        } catch (FileNotFoundException e) {
            System.out.printf("Error: cannot read the knownGeneCanonical file %s", knowgenecanonical);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 8. OMIM mim_pheno.txt file
        try {
            File mimpheno = new File(mim_pheno);
            try (BufferedReader reader_mim_pheno = new BufferedReader(new FileReader(mimpheno))) {
                String line2;
                while ((line2 = reader_mim_pheno.readLine()) != null) {
                    String[] cls2 = line2.split(" ");
                    if (cls2.length > 1) {
                        key = cls2[0];
                        value = cls2[1];
                        mim_pheno_dict.put(key, value);

                    }
                }
            }

        } catch (IOException e) {
            System.out.printf("Error: cannot read the MIM file %s", mim_pheno);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 9. OMIM mim_orpha.txt file
        try {
            File mimorpha = new File(mim_orpha);
            try (BufferedReader reader_mim_orpha = new BufferedReader(new FileReader(mimorpha))) {
                String line2;
                while ((line2 = reader_mim_orpha.readLine()) != null) {
                    String[] cls2 = line2.split(" ");
                    if (cls2.length > 1) {
                        key = cls2[0];
                        value = cls2[1];
                        mim_orpha_dict.put(key, value);
                    }
                }
            }

        } catch (IOException e) {
            System.out.printf("Error: cannot read the MIM file %s \n", mim_orpha);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 10. orpha.txt
        try {
            File Orpha = new File(orpha);
            try (BufferedReader reader_orpha = new BufferedReader(new FileReader(Orpha))) {
                String line2;
                while ((line2 = reader_orpha.readLine()) != null) {
                    String[] cls2 = line2.split("\t");
                    key = cls2[0];
                    value = cls2[1];
                    orpha_dict.put(key, value);
                }
            }

        } catch (IOException e) {
            System.out.printf("Error: can\\'t read the Orpha file %s", orpha);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 11. cancer_pathways=%(database_vic)s/cancer_pathways.list
        try {
            File Cancerpathway = new File(cancer_pathways);
            try (BufferedReader reader_cancer_pathways = new BufferedReader(new FileReader(Cancerpathway))) {
                String line2;
                while ((line2 = reader_cancer_pathways.readLine()) != null) {
                    String[] cls2 = line2.split("\t");
                    key = cls2[0];
                    value = "1";
                    cancer_pathways_dict.put(key, value);
                }
            }

        } catch (IOException e) {
            System.out.printf("Error: can\\'t read the cancer_pathways genes file %s", cancer_pathways);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 12. cancer_genes.list
        try {
            File cancergenes = new File(cancer_genes);
            try (BufferedReader reader_cancer_genes = new BufferedReader(new FileReader(cancergenes))) {
                String line2;
                while ((line2 = reader_cancer_genes.readLine()) != null) {
                    String[] cls2 = line2.split("\t");
                    key = cls2[0];
                    value = "1";
                    cancer_genes_dict.put(key, value);

                }
            }

        } catch (IOException e) {
            System.out.printf("Error: can\\'t read the cancer gene file %s \n", cancer_genes);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }

        // 13. VIC.cancer.types
        try {
            File cancertype = new File(cancer_types);
            try (BufferedReader reader_cancertype = new BufferedReader(new FileReader(cancertype))) {
                String line2;
                while ((line2 = reader_cancertype.readLine()) != null) {
                    String[] cls2 = line2.split("/");
                    key = cls2[1];
                    value = cls2[0];
// System.out.println("key and value of cancer types: " + key + "==" + value);
                    cancer_types_dict.put(key, value);
                }
            }

        } catch (IOException e) {
            System.out.printf("Error: can\\'t read the cancer_pathways genes file %s \n", cancer_types);
            System.out.println("Error: Please download it from the source website");
            System.exit(0);
        }
    }

    public static void Check_Downdb() throws IOException, InterruptedException {
        String cmd, dataset_file;
        File humandbpath = new File(database_humandb);
        if (!humandbpath.exists()) {
            humandbpath.mkdirs();
        } else {
            System.out.printf("INFO: The folder of %s is already created!\n", humandbpath);
        }

        // MY: updated to dbnsfp35a, clinvar_20190305, gnomad211_exome
        String ds = "refGene esp6500siv2_all 1000g2015aug avsnp150 " + dbnsfpVersion + " " + clinvarVersion + " exac03 dbscsnv11 dbnsfp31a_interpro ensGene knownGene " + cosmicVersion + " icgc21 " + gnomadVersion + "";
        for (String dbs : ds.split(" ")) {
            String file_name = dbs;
            if (dbs.equals("1000g2015aug")) {
                file_name = "ALL.sites.2015_08";
            }
            dataset_file = database_humandb + "/" + buildver + "_" + file_name + ".txt";
            File datasetfile = new File(dataset_file);
            if (!dbs.equals("rmsk")) {
                cmd = "perl " + annotatevariation + " -buildver " + buildver + " -downdb -webfrom annovar " + file_name + " " + database_humandb;

            } else {
                cmd = "perl " + annotatevariation + " -buildver " + buildver + " -downdb " + file_name + " " + database_humandb;

            }

            if (!datasetfile.isFile()) {
                if (dbs.equals("1000g2015aug")) {
                    file_name = "1000g2015aug";
                    dataset_file = database_humandb + "/" + buildver + "_" + file_name + ".txt";
                    cmd = "perl " + annotatevariation + "-buildver " + buildver + "-downdb -webfrom annovar " + file_name + " " + database_humandb;

                }
                System.out.printf("Warning: The Annovar dataset file of %s is not in %s, begin to download this %s ...\n", dbs, database_humandb, dataset_file);
                System.out.println(cmd);
                Process process = Runtime.getRuntime().exec(cmd);
                process.waitFor();
            }
        }
    }

    public static void Check_Input() throws IOException, InterruptedException {
        String cmd;
        if (inputtype.toLowerCase().equals("avinput")) {
            return;
        }

        File convert2annovar = new File(convert_annovar);
        if (inputtype.toLowerCase().equals("vcf")) {
            if (convert2annovar.isFile()) {
                cmd = "perl " + convert2annovar + " -format vcf4 " + inputfile + " -outfile " + inputfile + ".avinput"; //"-filter pass"@huyue REMOVED BY MY
                System.out.printf("INFO: Begin to convert your vcf file of %s to AVinput of Annovar ...\n", inputfile);
                System.out.println(cmd);
                Process process = Runtime.getRuntime().exec(cmd);
                InputStream ins = process.getInputStream();
                System.getProperty(cmd);
                process.waitFor();

            } else {
                System.out.printf("Error: The Annovar file [ %s ] is not here, please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar \n ", convert_annovar);
                System.exit(0);
            }
        }

        if (inputtype.toLowerCase().equals("vcf_m")) {
            if (convert2annovar.isFile()) {
                cmd = "perl " + convert2annovar + " -format vcf4 " + inputfile + " --allsample --outfile " + outputfile;
                System.out.printf("INFO: Begin to convert your vcf file with multiple samples of %s to AVinput of Annovar with All.raw.highqc.vcf.<samplename>.avinput...\n", inputfile);
                System.out.println("Warning: Please pay attention that the sample names in VCF file should contain letters/numners only, otherwise the converting may be failure!");
                System.out.println(cmd);
                Process process = Runtime.getRuntime().exec(cmd);
                System.getProperty(cmd);
                process.waitFor();
            } else {
                System.out.printf("Error: The Annovar file [ %s ] is not here, please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar", convert_annovar);
                System.exit(0);
            }
        }
    }

    public static void Check_annovar_result() throws IOException, InterruptedException {
        String annovar_options = "--otherinfo";
        String cmd;
        File table_annovar = new File(tableannovar);
        if (!table_annovar.isFile()) {
            System.out.printf("Error: The Annovar file [ %s ] is not here, please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar", tableannovar);
            System.exit(0);
        }

        switch (inputtype.toLowerCase()) {
            case "avinput": {
                cmd = "perl " + tableannovar + " " + inputfile + " " + database_humandb + " -buildver " + buildver + " -outfile " + outputfile
                        + " -remove -protocol refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp150," + dbnsfpVersion + ",dbscsnv11,"
                        + "dbnsfp31a_interpro," + clinvarVersion + "," + cosmicVersion + ",icgc21," + gnomadVersion + " -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f -nastring . " + annovar_options;
                System.out.println(cmd);
                Process process = Runtime.getRuntime().exec(cmd);
                System.getProperty(cmd);
                process.waitFor();
                break;
            }
            case "vcf": {
                cmd = "perl " + tableannovar + " " + inputfile + ".avinput " + database_humandb + " -buildver " + buildver + " -remove -out " + outputfile
                        + " -protocol refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp150," + dbnsfpVersion + ",dbscsnv11,"
                        + "dbnsfp31a_interpro," + clinvarVersion + "," + cosmicVersion + ",icgc21," + gnomadVersion + " -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f -nastring . " + annovar_options;
                System.out.println(cmd);
                Process process = Runtime.getRuntime().exec(cmd);
                System.getProperty(cmd);
                process.waitFor();
                break;
            }
            case "vcf_m":
                int o1 = -1;
                if (outputfile.contains("/")) {
                    o1 = outputfile.lastIndexOf("/");
                }
                String oname = outputfile.substring(o1 + 1);
                String opath = outputfile.substring(0, o1 + 1);
                if (opath.equals("")) {
                    opath = "./";
                }

                File output_path = new File(opath);
// System.out.println(oname + " for vcf_m");
// File outfile = new File(outputfile);
// File outputfile_path = new File(outfile.getAbsoluteFile().getParent());
                String[] f = output_path.list((File dir, String name) -> name.endsWith(".avinput") && name.startsWith(oname));
// System.out.println("=======" + f[0] + "======");
                for (String f1 : f) {
                    System.out.printf("INFO: Begin to annotate sample file of %s ....", f1);
                    String new_outfile = f1.replaceAll(".avinput", "");
                    cmd = "perl " + tableannovar + " " + opath + f1 + " " + database_humandb + " -buildver " + buildver + " -out " + opath + new_outfile
                            + " -remove -protocol refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp150," + dbnsfpVersion + ",dbscsnv11,"
                            + "dbnsfp31a_interpro," + clinvarVersion + "," + cosmicVersion + ",icgc21," + gnomadVersion + " -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f -nastring . " + annovar_options;
                    System.out.println(cmd);
                    Process process = Runtime.getRuntime().exec(cmd);
                    System.getProperty(cmd);
                    process.waitFor();
                }
                break;
            default:
                break;
        }
    }

    // def check_genes(anvfile):
    public static int Check_genes(String _anvfile) throws FileNotFoundException, IOException {
        int sum = 0;
        String grlfile = _anvfile + ".grl_p";
        File anvfile = new File(_anvfile);

        try (BufferedReader fh = new BufferedReader(new FileReader(anvfile))) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(grlfile))) {
                String line;
                int otherinfo_pos = 1;
                String gene_name;
                String line_out;
                while ((line = fh.readLine()) != null) {
                    String[] cls = line.split("\t");
                    if (cls.length > 1) {
                        if (sum == 0 && paras.get("otherinfo").toString().toLowerCase().contains("true")) {
                            for (int i = 0; i < cls.length; i++) {
                                if (cls[i].toLowerCase().contains("otherinfo")) {
                                    otherinfo_pos = i;
                                }
                            } // for
                        } // otherinfo = true
                        gene_name = cls[6];
                        if (cls[6].equals("Gene.refGene")) {
                            gene_name = "Gene";
                        }

                        sum = sum + 1;
                        for (String gg : gene_name.split(";")) {
                            if (!paras.get("otherinfo").toString().toLowerCase().contains("true")) {
                                line_out = line + "\t" + gg;
                            } else {
                                line_out = cls[0];
                                for (int i = 1; i < cls.length; i++) {
                                    if (cls.length > otherinfo_pos) {
                                        if (i != otherinfo_pos) {
                                            line_out = line_out + "\t" + cls[i];
                                        } else if (i == otherinfo_pos) {
                                            line_out = line_out + "\t" + gg + "\t" + cls[i];
                                        }
                                    } else {
                                        int o = otherinfo_pos - 1;
                                        if (i != o) {
                                            line_out = line_out + "\t" + cls[i];
                                        } else if (i == o) {
                                            line_out = line_out + "\t" + cls[i] + "\t" + gg;
                                        }
                                    } // else
                                } // for
                                if (sum > 1) {
                                    line_out = line_out.replace("(?i)chr", "");
                                }

                            } // else
                            writer.write(line_out + "\t\n");
                        } // for
                    } // cls.length
                } // while
            }
            fh.close();
        } catch (IOException e) {
            e.getStackTrace();
            System.out.printf("Error: can\'t read/write the annovar output file %s %s \n", anvfile, grlfile);
            System.exit(0);
        }

        return sum;

    }

    // def classfy(CBP, Allels_flgs, cls):
    public static String ClassFy(int[] CBP, Map Alleles_flgs, String[] cls) {

        String[] BPS = {"Strong clinical significance", "Potential clinical significance", "Benign/Likely benign", "Uncertain significance"};
        int BPS_out = 3;

// if (CBP[0] == 2 && CBP[3] == 1 && CBP[6] == 1 && CBP[7] == 2 && CBP[8] == 2 && CBP[9] == 2 && CBP[10] == 2) {
// BPS_out = 0;
// } else if (CBP[0] > 0 && CBP[3] > 0 && CBP[6] > 0 && CBP[7] == 2 && CBP[8] > 0 && CBP[9] == 2 && CBP[10] > 0) {
// BPS_out = 1;
// } else if (CBP[0] == 0 && CBP[3] == 0 && CBP[6] == 0 && CBP[7] == 0 && CBP[8] == 0 && CBP[9] == 0 && CBP[10] == 0) {
// BPS_out = 3;
// } else if (CBP[0] == 0 && CBP[3] == 0 && CBP[6] == 0 && CBP[7] >= 0 && CBP[8] >= 0 && CBP[9] <= 1 && CBP[10] == 0) {
// BPS_out = 2;
// }
//
// if (CBP[0] == 2 && CBP[1] == 1 && CBP[4] == 1 && CBP[5] == 2 && CBP[6] == 2 && CBP[7] == 2 && CBP[8] == 2) {
// BPS_out = 0;
// } else if (CBP[0] > 0 && CBP[1] > 0 && CBP[4] > 0 && CBP[5] == 2 && CBP[6] > 0 && CBP[7] == 2 && CBP[8] > 0) {
// BPS_out = 1;
// } else if (CBP[0] == 0 && CBP[1] == 0 && CBP[4] == 0 && CBP[5] == 0 && CBP[6] == 0 && CBP[7] == 1 && CBP[8] == 0) {
// BPS_out = 3;
// } else if (CBP[0] == 0 && CBP[1] == 0 && CBP[4] == 0 && CBP[5] <= 1 && CBP[6] == 0 && CBP[7] == 0 && CBP[8] == 0) {
// BPS_out = 2;
// }
//
// // BY MY 
// int CBP_57 = CBP[5] + CBP[7];
// if (CBP[0] == 2 && CBP[1] == 1 && CBP[4] == 1 && CBP[6] == 2 && CBP[8] == 2 && CBP_57 > 3) {
// BPS_out = 0;//strong evidence
// } else if (CBP[0] > 0 && CBP[1] > 0 && CBP[4] > 0 && CBP[6] > 0 && CBP[8] > 0 && CBP_57 > 1) {
// if (CBP[5] != 1 || CBP[7] != 0) {
// BPS_out = 1;//potential
// }
// } else if (CBP[0] == 0 && CBP[1] == 0 && CBP[4] == 0 && CBP[5] <= 1 && CBP[6] == 0 && CBP[7] <= 1 && CBP[8] == 0) {
// BPS_out = 2;//B/LB
// } else {
// BPS_out = 3;//VUS
// }
// 
        // BY HY
        int CBP_57 = CBP[5] + CBP[7];
        if (CBP[0] == 2 && CBP[1] == 1 && CBP[4] == 1 && CBP[6] == 2 && CBP[8] == 2 && CBP_57 > 3) {
            BPS_out = 0;
        } else if (CBP[0] > 0 && CBP[1] > 0 && CBP[4] > 0 && CBP[6] > 0 && CBP[8] > 0 && CBP_57 > 1) {
            if (CBP[5] != 1 || CBP[7] != 1) {
                BPS_out = 1;
            }
        } else if (CBP[0] == 0 && CBP[1] == 0 && CBP[4] == 0 && CBP[5] <= 1 && CBP[6] == 0 && CBP[7] <= 1 && CBP[8] == 0) {
            BPS_out = 2;
        } else {
            BPS_out = 3;
        }

        return BPS[BPS_out];
    }

    public static int Check_Thera(String line, Map Funcanno_flgs) throws FileNotFoundException, IOException {
        int level = 0;
        String[] cls = line.split("\t");

        int cl = Integer.parseInt(Funcanno_flgs.get("Otherinfo").toString());

        if (cl < cls.length && cl != 0) {
            String clss = cls[cl];
            String[] clstt = clss.split(";");

            if (paras.containsKey("cancer_type")) {
                cancer_type = paras.get("cancer_type").toString();//assign cancer_type by using '-cancer_type option', e.g: AA (always use short cuts)
            } else if (!paras.containsKey("cancer_type")) {//assign cancer_types by putting cancer types in the otherinfo column in the Avinput file, use ";" as seperator, VIC will take the first cancer type in the list.
                if (clstt[0].length() > 0) {
                    cancer_type = clstt[0];
                }
            }
        } else {
            cancer_type = "CANCER";
        }

        cancer_type = cancer_type.toUpperCase();
        int getr = Integer.parseInt(Funcanno_flgs.get("Gene").toString());
        String gene_tr = cls[getr];

        int fr = Integer.parseInt(Funcanno_flgs.get("Func.refGene").toString());
        String func = cls[fr];
        int erg = Integer.parseInt(Funcanno_flgs.get("ExonicFunc.refGene").toString());
        String exon_func = cls[erg];

        String line_tmp = gene_tr + " " + func + " " + exon_func + " " + cancer_type;

        int aarf = Integer.parseInt(Funcanno_flgs.get("AAChange.refGene").toString());
        String line_tmp2 = cls[aarf];

        String[] line_tmp2_sp = line_tmp2.split(", ");
        String marker_key0; // = gene_tr + "_" + "mut_any";

        String level_AB = "0";
        String level_CD = "0";
        String level_cancer = "0";
        String level_CD_cancer = "0";
        for (String cls0 : line_tmp2_sp) {
            String[] cls0_1 = cls0.split(":");
            if (cls0.length() > 1 && line_tmp2_sp.length > 0 && cls0_1.length == 5) {
// cls0_1 = cls0.split(":");
                String gene = cls0_1[0];
                String transc = cls0_1[1];
                String exon = cls0_1[2].replaceAll("exon", "");
                String cdna_change = cls0_1[3].replaceAll("c.", "");//c.G38A
                String amino_change = cls0_1[4].replaceAll("p.", "");//p.G13D

                int ltt = amino_change.length();
                String codon_num = amino_change.substring(1, ltt - 1);

                String marker_key = gene + "_" + amino_change;
                marker_key0 = gene + "_mut_any";

                String marker_key1 = gene + "_exon_" + exon + "_any";
                String marker_key2 = gene + "_" + "codon_" + codon_num + "_any";
                String marker_key00 = "";
                String marker_key11 = "";
                String marker_key22 = "";

                if (exon_func.equals("nonsynonymous SNV")) {
                    marker_key00 = gene + "_missense_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_missense";
                    marker_key22 = gene + "_codon_" + codon_num + "_missense";

                } else if (exon_func.contains("frameshift") && !exon_func.contains("nonframe")) {
                    marker_key00 = gene + "_frameshift_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_frameshift";
                    marker_key22 = gene + "_codon_" + codon_num + "_frameshift";

                } else if (exon_func.contains("stopgain") || exon_func.contains("stoploss")) {
                    marker_key00 = gene + "_nonsense_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_nonsense";
                    marker_key22 = gene + "_codon_" + codon_num + "_nonsense";

                } else if (exon_func.contains("deletion")) {
                    marker_key00 = gene + "_deletion_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_deletion";
                    marker_key22 = gene + "_codon_" + codon_num + "_deletion";
                } else if (exon_func.contains("insertion")) {
                    marker_key00 = gene + "_insertion_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_insertion";
                    marker_key22 = gene + "_codon_" + codon_num + "_insertion";
                }

                String add_list = (String) add_markers_dict.get(marker_key) + add_markers_dict.get(marker_key0) + add_markers_dict.get(marker_key00) + add_markers_dict.get(marker_key1)
                        + add_markers_dict.get(marker_key11) + add_markers_dict.get(marker_key2) + add_markers_dict.get(marker_key22);

                String cgi_list = (String) cgi_markers_dict.get(marker_key) + cgi_markers_dict.get(marker_key0) + cgi_markers_dict.get(marker_key00) + cgi_markers_dict.get(marker_key1)
                        + cgi_markers_dict.get(marker_key11) + cgi_markers_dict.get(marker_key2) + cgi_markers_dict.get(marker_key22);
// System.out.println("=====cgi_list is: " + cgi_list);
// System.out.println("=====###$$%add_list is " + add_list);
// System.out.println("=====###$$%cgi_list is " + cgi_list);
                String level_AB_tmp; // = "0";
                String level_CD_tmp; // = "0";
                int level_tmp; // = 0;
                for (String i : add_list.split(", ")) {
                    if (i.length() <= 8 && i.length() > 0) {
                        level_CD_tmp = "1";
                        level_tmp = 1;
                        if (Integer.valueOf(level_CD_tmp) > Integer.valueOf(level_CD)) {
                            level_CD = level_CD_tmp;
                        }
                        if (level_tmp > level) {
                            level = level_tmp;
                        }
                    }
                }

                try (InputStream in = new FileInputStream(cgi_markers)) {
                    byte[] b = new byte[in.available()];
                    in.read(b);
                    String[] cgiline = new String(b).split("\n");
                    for (String i : cgi_list.split(", ")) {
                        if (0 < i.length() && i.length() <= 8) {
                            int pos = Integer.parseInt(i);
                            String cgid = cgiline[pos];

                            String[] cgidd = cgid.split("\t");
// System.out.println("=====cgidd is: " + Arrays.toString(cgidd));

                            if (cgidd[9].toLowerCase().contains("pre-clinical") || cgidd[9].toLowerCase().contains("case report")) {
                                level_tmp = 1;
                                level_CD_tmp = "1";
                                if (Integer.valueOf(level_CD_tmp) > Integer.valueOf(level_CD)) {
                                    level_CD = level_CD_tmp;
                                }
                                if (level_tmp > level) {
                                    level = level_tmp;
                                }
                            }
                            if (cgidd[9].toLowerCase().contains("trial") || cgidd[5].toLowerCase().contains("trial")) {
                                level_tmp = 1;
                                level_CD_tmp = "1";
                                if (Integer.valueOf(level_CD_tmp) > Integer.valueOf(level_CD)) {
                                    level_CD = level_CD_tmp;
                                }
                                if (level_tmp > level) {
                                    level = level_tmp;
                                }
                            }
                            if (cgidd[9].contains("guidelines") || cgidd[5].contains("approved")) {
                                level_tmp = 2;
                                level_AB_tmp = "1";
                                if (Integer.valueOf(level_AB_tmp) > Integer.valueOf(level_AB)) {
                                    level_AB = level_AB_tmp;
                                }
                                if (level_tmp > level) {
                                    level = level_tmp;
                                }
                            }
                            if (!cgidd[14].toUpperCase().contains(cancer_type)) {
                                level_cancer = "1";
                                level_CD_cancer = "1";
                            } else if (cgidd[14].toUpperCase().contains(cancer_type)) {
                                level_cancer = "2";
                                level_CD_cancer = "2";
                            }

                            if (Integer.valueOf(level_cancer) > level) {
                                level = Integer.valueOf(level_cancer);
                            }
                        }
                    }
                    if (level_AB.equals("1")) {
                        level = 2;
                    }
                } catch (Exception e) {
                    System.out.println(e);
                }
            }
        }

        return level;
    }

    public static String Check_Thera_out(String line, Map Funcanno_flgs) throws FileNotFoundException, IOException {
        String thera_out = "";
        String[] cls = line.split("\t");

        int cl = Integer.parseInt(Funcanno_flgs.get("Otherinfo").toString());

        if (cl < cls.length) {
            String clss = cls[cl];
            String[] clstt = clss.split(";");

            if (paras.containsKey("cancer_type")) {
                cancer_type = paras.get("cancer_type").toString();
            } else if (!paras.containsKey("cancer_type")) {
                if (clstt[0].length() > 0) {
                    cancer_type = clstt[0];
                }
            }
        } else {
            cancer_type = "CANCER";
        }

        int getr = Integer.parseInt(Funcanno_flgs.get("Gene").toString());
        String gene_tr = cls[getr];

        int fr = Integer.parseInt(Funcanno_flgs.get("Func.refGene").toString());
        String func = cls[fr];
        int erg = Integer.parseInt(Funcanno_flgs.get("ExonicFunc.refGene").toString());
        String exon_func = cls[erg];

// String line_tmp = gene_tr + " " + func + " " + exon_func + " " + cancer_type;
        int aarf = Integer.parseInt(Funcanno_flgs.get("AAChange.refGene").toString());
        String line_tmp2 = cls[aarf];

        String[] line_tmp2_sp = line_tmp2.split(", ");
        String marker_key0; // = gene_tr + "_" + "mut_any";

        for (String cls0 : line_tmp2_sp) {
            String[] cls0_1 = cls0.split(":");
            if (cls0.length() > 1 && line_tmp2_sp.length > 0 && cls0_1.length == 5) {
                cls0_1 = cls0.split(":");
                String gene = cls0_1[0];
                String transc = cls0_1[1];
                String exon = cls0_1[2].replaceAll("exon", "");
                String cdna_change = cls0_1[3].replaceAll("c.", "");//c.G38A
                String amino_change = cls0_1[4].replaceAll("p.", "");//p.G13D

                int ltt = amino_change.length();
                String codon_num = amino_change.substring(1, ltt - 1);

                String marker_key = gene + "_" + amino_change;
                marker_key0 = gene + "_mut_any";

                String marker_key1 = gene + "_exon_" + exon + "_any";
                String marker_key2 = gene + "_" + "codon_" + codon_num + "_any";
                String marker_key00 = "";
                String marker_key11 = "";
                String marker_key22 = "";

                if (exon_func.equals("nonsynonymous SNV")) {
                    marker_key00 = gene + "_missense_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_missense";
                    marker_key22 = gene + "_codon_" + codon_num + "_missense";

                } else if (exon_func.contains("frameshift") && !exon_func.contains("nonframe")) {
                    marker_key00 = gene + "_frameshift_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_frameshift";
                    marker_key22 = gene + "_codon_" + codon_num + "_frameshift";

                } else if (exon_func.contains("stopgain") || exon_func.contains("stoploss")) {
                    marker_key00 = gene + "_nonsense_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_nonsense";
                    marker_key22 = gene + "_codon_" + codon_num + "_nonsense";

                } else if (exon_func.contains("deletion")) {
                    marker_key00 = gene + "_deletion_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_deletion";
                    marker_key22 = gene + "_codon_" + codon_num + "_deletion";
                } else if (exon_func.contains("insertion")) {
                    marker_key00 = gene + "_insertion_any";
                    marker_key11 = gene + "_" + "exon_" + exon + "_insertion";
                    marker_key22 = gene + "_codon_" + codon_num + "_insertion";
                }

                String add_list = (String) add_markers_dict.get(marker_key) + add_markers_dict.get(marker_key0) + add_markers_dict.get(marker_key00) + add_markers_dict.get(marker_key1)
                        + add_markers_dict.get(marker_key11) + add_markers_dict.get(marker_key2) + add_markers_dict.get(marker_key22);

                String cgi_list = (String) cgi_markers_dict.get(marker_key) + cgi_markers_dict.get(marker_key0) + cgi_markers_dict.get(marker_key00) + cgi_markers_dict.get(marker_key1)
                        + cgi_markers_dict.get(marker_key11) + cgi_markers_dict.get(marker_key2) + cgi_markers_dict.get(marker_key22);
// System.out.println("=====cgi_list is: " + cgi_list);
// String level_AB = "0";
// String level_CD = "0";

                try (InputStream in = new FileInputStream(cgi_markers)) {
                    byte[] b = new byte[in.available()];
                    in.read(b);
                    String[] cgiline = new String(b).split("\n");
                    for (String i : cgi_list.split(", ")) {
                        if (0 < i.length() && i.length() <= 8) {
                            int pos = Integer.parseInt(i);
                            String cgid = cgiline[pos];

                            String[] cgidd = cgid.split("\t");
// System.out.println("=====cgidd is: " + Arrays.toString(cgidd));
                            thera_out = Arrays.toString(cgidd);
                        }
                    }
                } catch (Exception e) {
                    System.out.println(e);
                }
            }
        }

        return thera_out;
    }

    public static int Check_Mut(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        /* Mutation type:
 1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
 1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
 0 Functionally unknown; mostly missense, in-frame indels; less commonly, other types
 0 Functionally benign or unknown; mostly missense; less commonly, other types
         */

        String[] cls = line.split("\t");
        List<String> funcs_tmp = new ArrayList<>(asList("nonsynonymous", "missense", "nonsense", "frameshift", "splic", "stopgain", "stoplost", "CNV", "fusion"));
        int rg = Integer.parseInt(Funcanno_flgs.get("Func.refGene").toString());
        String func = cls[rg];
        int erg = Integer.parseInt(Funcanno_flgs.get("ExonicFunc.refGene").toString());
        String exon_func = cls[erg];
        int getr = Integer.parseInt(Funcanno_flgs.get("Gene").toString());

        //String funcs_tmp2 = "nonframe";
        String line_tmp = func + " " + exon_func;

        int Mut = 0;
        int VS_t1 = 0;
        int VS_t2 = 0;
        int VS_t3 = 0;
        double dbscSNV_cutoff = 0.6;

        for (String fc : funcs_tmp) {
            // if (line_tmp.contains(fc) && !line_tmp.contains(funcs_tmp2)) {
            if (line_tmp.contains(fc)) {
                VS_t1 = 1;
                break;
            }
        }

        if (lof_genes_dict.containsKey(cls[getr])) {
            if (lof_genes_dict.get(cls[getr]).equals("1")) {
                VS_t2 = 1;
            }
        } else {
            VS_t2 = 0;
        }

        int rf = Integer.parseInt(Funcanno_flgs.get("dbscSNV_RF_SCORE").toString());
        int ada = Integer.parseInt(Funcanno_flgs.get("dbscSNV_ADA_SCORE").toString());

        if (!cls[rf].equals(".") && !cls[ada].equals(".")) {
            double rf_score = Double.parseDouble(cls[rf]);
            double ada_score = Double.parseDouble(cls[ada]);
            if (rf_score > dbscSNV_cutoff || ada_score > dbscSNV_cutoff) {
                VS_t3 = 1;
            } else {
                VS_t2 = 0;
            }
        }

        if (VS_t1 != 0 && VS_t2 != 0) {
            Mut = 1;
        }
        if (VS_t3 != 0 && VS_t2 != 0) {
            Mut = 1;
        }

        return Mut;
    }

    public static int Check_VF(String line, Map Funcanno_flgs, Map Alleles_flgs, Map lof_genes_dict) {
        /* Variant frequencies
 1 Mostly mosaic
 1 Mostly mosaic
 0 Mosaic or nonmosaic
 0 Mostly nonmosaic (VAF, approximately 50% or 100%)
         */
        int VF = -1;
        return VF;
    }

    public static int Check_PotG(String line, Map Funcanno_flgs, Map Alleles_flgs, Map lof_genes_dict) {
        /* Potential germline
 Mostly nonmosaic (VAF approximately 50% or 100%)
 Mostly nonmosaic (VAF approximately 50% or 100%)
 Mostly nonmosaic (VAF approximately 50% or 100%)
 Mostly nonmosaic (VAF, approximately 50% or 100%)
         */
        int PotG = -1;

        return PotG;
    }

    public static int Check_PopD(String line, Map Freqs_flgs, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        /* Population database: ESP, dbSNP, 1000Genome, ExAC
 1 Absent or extremely low MAF
 1 Absent or extremely low MAF
 1 Absent or extremely low MAF
 0 MAF>1% in the general population; or high MAF in some ethnic populations
         */

        double MAF_cutoff = 0.005;
        int PopD = 0;

        String[] cls = line.split("\t");
        String[] Freqs_4pops = {"1000g2015aug_all", "esp6500siv2_all", "ExAC_ALL", "AF"};//AF in gnomAD
        int tt = 1;

        for (String kkey : Freqs_4pops) {
            int k = Integer.parseInt(Freqs_flgs.get(kkey).toString());
            if (!cls[k].equals(".")) {
                tt = tt * 0;
            }
        }

        if (tt == 1) {
            PopD = 1;
        }

        for (String popkey : Freqs_4pops) {
            try {
                int popvalue = Integer.parseInt(Freqs_flgs.get(popkey).toString());
                if (!cls[popvalue].equals(".")) {
                    if (Double.parseDouble(cls[popvalue]) > 0.01) {
                        PopD = 0;
                    } else if (Double.parseDouble(cls[popvalue]) < MAF_cutoff) {
                        PopD = 1;
                    }
                }

            } catch (ValueException e) {
            }

        }
        return PopD;
    }

// 
    public static int Check_GermD(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        int GermD = 0;
        String[] cls = line.split("\t");
        int cli = Integer.parseInt(Funcanno_flgs.get("CLNSIG").toString());
        String line_tmp2 = cls[cli];
        if (line_tmp2.contains("onflicting")) {
            GermD = 0;
        } else if (!line_tmp2.contains("enign") && line_tmp2.contains("athogenic")) {
            GermD = 2;
        } else if (!line_tmp2.contains("enign") && line_tmp2.contains("drug")) {
            GermD = 2;
        } else if (line_tmp2.contains("ikely_benign") || line_tmp2.contains("enign")) {
            GermD = 1;
        } else if (line_tmp2.contains("enign") && line_tmp2.contains("athogenic")) {
            GermD = 0;
        } else if (line_tmp2.contains("ncertain")) {
            GermD = 0;
        }
        return GermD;
    }

    //
    public static int Check_SomD(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        /*
 Somatic database: COSMIC, My Cancer Genome, TCGA
 2 Most likely present
 1 Likely present
 0 Absent or present without association to specific tumors (potential germline VUS); present but in very few cases
 0 Absent or present without association to specific tumors (potential rare germline polymorphism)
         */
        int SomD = 0;
        String[] cls = line.split("\t");
        int cos = Integer.parseInt(Funcanno_flgs.get(cosmicVersion).toString());
        int icgc = Integer.parseInt(Funcanno_flgs.get("ICGC_Id").toString());

        if (!cls[cos].equals(".") || !cls[icgc].equals(".")) {
            SomD = 1;
        }
        if (cls[cos].equals(".") && cls[icgc].equals(".")) {
            SomD = 0;
        }
        if (!cls[cos].equals(".") && !cls[icgc].equals(".")) {
            SomD = 2;
        }
        return SomD;
    }

    // def check_PreP(line, Funcanno_flgs, Allels_flgs, lof_genes_dict)
    public static int Check_PreP(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        /*
 Predictive software: SIFT, PolyPhen2, MutTaster, CADD, MetaSVM
 2 Mostly damaging; information to be used for reference only
 0 Variable; information to be used for reference only
 1 Mostly benign; information to be used for reference only
         */
        double sift_cutoff = 0.05;
        int metasvm_cutoff = 0;
        int cutoff_conserv = 2;

        int dam = 0;
        int var = 0;
        int ben = 0;
        int PreP = 0;

        String[] cls = line.split("\t");
        int svm = Integer.parseInt(Funcanno_flgs.get("MetaSVM_score").toString());
        int sift = Integer.parseInt(Funcanno_flgs.get("SIFT_score").toString());

        if (!cls[svm].equals(".")) {
            if (Double.parseDouble(cls[svm]) < metasvm_cutoff) {
                ben = ben + 1;
            } else {
                dam = dam + 1;
            }
        } else {
            var = var + 1;
        }

        if (!cls[sift].equals(".")) {
            if (Double.parseDouble(cls[sift]) >= sift_cutoff) {
                ben = ben + 1;
            } else {
                dam = dam + 1;
            }
        } else {
            var = var + 1;
        }

        int ph = Integer.parseInt(Funcanno_flgs.get("Polyphen2_HDIV_pred").toString());
        switch (cls[ph]) {
            case "P":
            case "D":
                dam = dam + 1;
                break;
            case "B":
                ben = ben + 1;
                break;
            case ".":
                var = var + 1;
                break;
            default:
                break;
        }

        int lr = Integer.parseInt(Funcanno_flgs.get("MetaLR_pred").toString());
        switch (cls[lr]) {
            case "D":
                dam = dam + 1;
                break;
            case "T":
                ben = ben + 1;
                break;
            case ".":
                var = var + 1;
                break;
            default:
                break;
        }

        int fa = Integer.parseInt(Funcanno_flgs.get("FATHMM_pred").toString());
        switch (cls[fa]) {
            case "D":
                dam = dam + 1;
                break;
            case "T":
                ben = ben + 1;
                break;
            case ".":
                var = var + 1;
                break;
            default:
                break;
        }

        int taster = Integer.parseInt(Funcanno_flgs.get("MutationTaster_pred").toString());
        switch (cls[taster]) {
            case "A":
            case "D":
                dam = dam + 1;
                break;
            case "P":
                ben = ben + 1;
                break;
            case ".":
                var = var + 1;
                break;
            default:
                break;
        }

        int gerp = Integer.parseInt(Funcanno_flgs.get("GERP++_RS").toString());
        if (cls[gerp].equals(".")) {
            var = var + 1;
        } else {

            double gerp_score = Double.parseDouble(cls[gerp]);
            if (gerp_score >= cutoff_conserv) {
                dam = dam + 1;
            } else {
                ben = ben + 1;
            }
        }

        if (dam > 3) {
            PreP = 2;
        }
        if (ben > 3) {
            PreP = 0;
        }
        if (var > 3) {
            PreP = 1;
        }
        if (dam == ben) {
            PreP = 1;
        }
        return PreP;
    }

    public static int Check_Path(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        int Path = 0;
        String cls[] = line.split("\t");

        int getr = Integer.parseInt(Funcanno_flgs.get("Gene").toString());
        String gene_tr = cls[getr];

        if (cancer_pathways_dict.containsKey(gene_tr)) {
            if (cancer_pathways_dict.get(gene_tr) == "1") {
                Path = 1;
            }
        }

        if (cancer_genes_dict.containsKey(gene_tr)) {
            if (cancer_genes_dict.get(gene_tr) == "1") {
                Path = 2;
            }
        }

        return Path;
    }

    public static int Check_Pubs(String line, Map Funcanno_flgs, Map Allels_flgs, Map lof_genes_dict) {
        /*
 Publications: functional study, population study, other
 Therapeutic/Diagnostic/Prognostic: reported evidence with consensus
 Therapeutic: evidence of using FDA-approved therapies for different tumor types; phase 2 or 3 clinical trials for investigational therapies; Diagnostic/Prognostic: multiple lines of reported evidence without consensus
 None or no convincing evidence to determine clinical/biological significance
 Reported evidence supportive of benign/likely benign; or none
         */
        int Pubs = -1;
        return Pubs;
    }

    public static String Assign(String BP, String line, Map Freqs_flgs, Map Funcanno_flgs, Map Allels_flgs) throws FileNotFoundException, IOException {
        int[] CBP = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        String BP_out = null;
        int Therapeutic = Check_Thera(line, Funcanno_flgs);
        CBP[0] = Therapeutic;

        int Mutation_type = Check_Mut(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[1] = Mutation_type;

        int Variant_freq = Check_VF(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[2] = Variant_freq;

        int Potential_germ = Check_PotG(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[3] = Potential_germ;

        int Population_data = Check_PopD(line, Freqs_flgs, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[4] = Population_data;

        int Germline_data = Check_GermD(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[5] = Germline_data;

        int Somatic_data = Check_SomD(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[6] = Somatic_data;

        int Predict_pathogenic = Check_PreP(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[7] = Predict_pathogenic;

        int Pathway_invol = Check_Path(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[8] = Pathway_invol;

        int Publications = Check_Pubs(line, Funcanno_flgs, Allels_flgs, lof_genes_dict);
        CBP[9] = Publications;

        String[] cls = line.split("\t");

        if (!evidencefile.equals(" ")) {
            File evidence_file = new File(evidencefile);
            if (evidence_file.isFile()) {
                int chr = (int) Allels_flgs.get("Chr");
                int start = (int) Allels_flgs.get("Start");
                int ref = (int) Allels_flgs.get("Ref");
                int alt = (int) Allels_flgs.get("Alt");
                String keys = cls[chr] + "_" + cls[start] + "_" + cls[ref] + "_" + cls[alt];
// System.out.println("======target variant is " + keys);
                keys = keys.replaceAll("(?!)chr", "");
                if (user_evidence_dict.get(keys) != null) {
                    String evds = (String) user_evidence_dict.get(keys.toUpperCase());
// System.out.println("======evidence are: " + evds);
                    for (String evd : evds.split(";")) {
                        String[] evd_t = evd.split("=");
// System.out.println(evd_t[0] + " " + evd_t[1]);
                        switch (evd_t[0]) {
                            case "CBP0":
                                CBP[0] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP1":
                                CBP[1] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP2":
                                CBP[2] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP3":
                                CBP[3] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP4":
                                CBP[4] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP5":
                                CBP[5] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP6":
                                CBP[6] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP7":
                                CBP[7] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP8":
                                CBP[8] = Integer.valueOf(evd_t[1]);
                                break;
                            case "CBP9":
                                CBP[9] = Integer.valueOf(evd_t[1]);
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
        }

        if (cls.length > 1) {
            BP_out = ClassFy(CBP, Allels_flgs, cls);
            String line_t = BP_out + " EVS=" + Arrays.toString(CBP);
            BP_out = line_t.replaceAll("-1", "NONE");
        }

        return BP_out;
    }

    public static void Search_key_index(String line, Map dict) {
        String[] cls = line.split("\t");
        for (Object kkey : dict.keySet()) {
            for (int i = 1; i < cls.length + 1; i++) {
                int ii = i - 1;
                if (kkey.toString().equals(cls[ii])) {
                    dict.put(kkey, ii);
                    break;
                }
            }
        }
    }

    // def my_inter_var_can(annovar_outfile):
    public static int My_Cancer_var_can(String _annovarfile) throws FileNotFoundException, IOException {

        String grlfile = _annovarfile + ".grl_p";
        String vicfile = _annovarfile + ".vic";
        line_sum = 0;

        Freqs_flgs.put("1000g2015aug_all", 0);
        Freqs_flgs.put("esp6500siv2_all", 0);
        Freqs_flgs.put("ExAC_ALL", 0);
        Freqs_flgs.put("ExAC_AFR", 0);
        Freqs_flgs.put("ExAC_AMR", 0);
        Freqs_flgs.put("ExAC_EAS", 0);
        Freqs_flgs.put("ExAC_FIN", 0);
        Freqs_flgs.put("ExAC_NFE", 0);
        Freqs_flgs.put("ExAC_OTH", 0);
        Freqs_flgs.put("ExAC_SAS", 0);
        Freqs_flgs.put("AF", 0);
        Freqs_flgs.put("AF_afr", 0);
        Freqs_flgs.put("AF_sas", 0);
        Freqs_flgs.put("AF_amr", 0);
        Freqs_flgs.put("AF_eas", 0);
        Freqs_flgs.put("AF_nfe", 0);
        Freqs_flgs.put("AF_fin", 0);
        Freqs_flgs.put("AF_asj", 0);
        Freqs_flgs.put("AF_oth", 0);

        Funcanno_flgs.put("Func.refGene", 0);
        Funcanno_flgs.put("ExonicFunc.refGene", 0);
        Funcanno_flgs.put("AAChange.refGene", 0);
        Funcanno_flgs.put("Gene", 0);
        Funcanno_flgs.put("Gene damage prediction (all disease-causing genes)", 0);
        Funcanno_flgs.put("CLNDBN", 0);
        Funcanno_flgs.put("CLNACC", 0);
        Funcanno_flgs.put("CLNDSDB", 0);
        Funcanno_flgs.put("dbscSNV_ADA_SCORE", 0);
        Funcanno_flgs.put("dbscSNV_RF_SCORE", 0);
        Funcanno_flgs.put("GERP++_RS", 0);
        Funcanno_flgs.put("LoFtool_percentile", 0);
        Funcanno_flgs.put("Interpro_domain", 0);
        Funcanno_flgs.put("rmsk", 0);
        Funcanno_flgs.put("SIFT_score", 0);
        Funcanno_flgs.put("phyloP46way_placental", 0);
        Funcanno_flgs.put("Gene.ensGene", 0);
        Funcanno_flgs.put("CLNSIG", 0);

        Funcanno_flgs.put("CADD_phred", 0);
        Funcanno_flgs.put("avsnp150", 0);
        Funcanno_flgs.put("AAChange.ensGene", 0);
        Funcanno_flgs.put("AAChange.knownGene", 0);
        Funcanno_flgs.put("MetaSVM_score", 0);
        Funcanno_flgs.put(cosmicVersion, 0);
        Funcanno_flgs.put("ICGC_Id", 0);
        Funcanno_flgs.put("ICGC_Occurrence", 0);
        Funcanno_flgs.put("Otherinfo", 0);
        Funcanno_flgs.put("Polyphen2_HDIV_pred", 0);
        Funcanno_flgs.put("MetaLR_pred", 0);
        Funcanno_flgs.put("MutationTaster_pred", 0);
        Funcanno_flgs.put("FATHMM_pred", 0);

        Allels_flgs.put("Chr", 0);
        Allels_flgs.put("Start", 0);
        Allels_flgs.put("End", 0);
        Allels_flgs.put("Ref", 0);
        Allels_flgs.put("Alt", 0);

        File fh = new File(grlfile);
        File fw = new File(vicfile);
        BufferedReader reader = new BufferedReader(new FileReader(fh));
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fw))) {
            line_sum = 0;
            String strs;
            System.out.println("Notice:Begin the variants interpretation by VIC......");
            if (paras.get("otherinfo").toString().toLowerCase().contains("true")) {
                writer.write("Chr\t" + "Start\t" + "End\t" + "Ref\t" + "Alt\t" + "Ref.Gene\t" + "Func.refGene\t" + "ExonicFunc.refGene\t"
                        + "Gene.ensGene\t" + "avsnp150\t" + "clinvar: " + "Clinvar\t" + "VIC: " + "VIC and Evidence\t" + "AAChange.ensGene\t" + "AAChange.refGene\t" + "AAChange.knownGene\t"
                        + "Freq_ExAC_ALL\t" + "Freq_esp6500siv2_all\t" + "Freq_1000g2015aug_all\t" + "Freq_gnomAD_AF\t" + "FATHMM_pred\t" + "CADD_phred\t"
                        + "SIFT_score\t" + "GERP++_RS\t" + "Polyphen2_HDIV_pred\t" + "dbscSNV_ADA_SCORE\t" + "dbscSNV_RF_SCORE\t"
                        + "Interpro_domain\t" + "MetaSVM_score\t" + "MutationTaster_pred\t" + "Freq_ExAC_POPs\t" + "Freq_gnomAD_exome_POPs\t" + "OMIM\t" + "Phenotype_MIM\t"
                        + "OrphaNumber\t" + "Orpha\t" + "cgi_info\t" + "civic_info\t" + "Otherinfo\n");
            } else {
                writer.write("Chr\t" + "Start\t" + "End\t" + "Ref\t" + "Alt\t" + "Ref.Gene\t" + "Func.refGene\t" + "ExonicFunc.refGene\t"
                        + "Gene.ensGene\t" + "avsnp150\t" + "clinvar: " + "Clinvar\t" + "VIC: " + "VIC and Evidence\t" + "AAChange.ensGene\t" + "AAChange.refGene\t" + "AAChange.knownGene\t"
                        + "Freq_ExAC_ALL\t" + "Freq_esp6500siv2_all\t" + "Freq_1000g2015aug_all\t" + "Freq_gnomAD_AF\t" + "FATHMM_pred\t" + "CADD_phred\t"
                        + "SIFT_score\t" + "GERP++_RS\t" + "Polyphen2_HDIV_pred\t" + "dbscSNV_ADA_SCORE\t" + "dbscSNV_RF_SCORE\t"
                        + "Interpro_domain\t" + "MetaSVM_score\t" + "MutationTaster_pred\t" + "Freq_ExAC_POPs\t" + "Freq_gnomAD_exome_POPs\t" + "OMIM\t" + "Phenotype_MIM\t"
                        + "OrphaNumber\t" + "Orpha\t" + "cgi_info\t" + "civic_info\n");
            }
            writer.flush();

            while ((strs = reader.readLine()) != null) {
                for (String line : strs.split("\n")) {
// System.out.println("line in strs@1786 is: " + line);
                    String BP = "UNK";
                    String clinvar_bp; // = "UNK";
                    String[] cls = line.split("\t");
                    String each_variant = cls[0] + "#" + cls[1] + "#" + cls[2] + "#" + cls[3] + "#" + cls[4];
                    each_variant = each_variant.replaceAll("(?i)chr", "");
// System.out.println("variant is: " + each_variant);
                    String civicdd = "";
                    if (civ_markers_dict.containsKey(each_variant)) {
                        civicdd = civ_markers_dict.get(each_variant).toString();
                    }
                    if (civicdd.equals("")) {
                        civicdd = ".";
                    }
                    if (cls.length < 2) {
                        break;
                    }
                    if (line_sum == 0) {
                        Search_key_index(line, Freqs_flgs);
                        Search_key_index(line, Funcanno_flgs);
                        Search_key_index(line, Allels_flgs);
                    } else {
                        int intclinsig = Integer.parseInt(Funcanno_flgs.get("CLNSIG").toString());
                        String line_tmp2 = cls[intclinsig];
                        if (!line_tmp2.equals(".")) {
                            String[] cls3 = line_tmp2.split(";");
                            clinvar_bp = cls3[0];
                        } else {
                            clinvar_bp = ".";
                        }
// System.out.println("=========Gene is:" + cls[Funcanno_flgs.get("Gene")] + "===");

                        // 03142019 by MQY
                        String vic_bp = "";
                        if (knownlist.equals(" ")) {
                            if (clinvar_bp.equals(".") && cls[Funcanno_flgs.get("Gene")].equals(".") && cls[Funcanno_flgs.get("Func.refGene")].equals(".")) {
                                vic_bp = ".";
                            } else {
                                vic_bp = Assign(BP, line, Freqs_flgs, Funcanno_flgs, Allels_flgs);
                            }
                        } else {
                            for (String known_variant : known_list) {
                                if (known_variant.contains(each_variant)) {
// System.out.println("==" + known_variant + "==");
                                    vic_bp = "User-specified as " + known_variant.split("#")[5];
                                    break;
                                }
                                if (clinvar_bp.equals(".") && cls[Funcanno_flgs.get("Gene")].equals(".") && cls[Funcanno_flgs.get("Func.refGene")].equals(".")) {
                                    vic_bp = ".";
                                } else {
                                    vic_bp = Assign(BP, line, Freqs_flgs, Funcanno_flgs, Allels_flgs);
                                }
                            }
                        }

                        String Freq_ExAC_POPs = "AFR:" + cls[Freqs_flgs.get("ExAC_AFR")] + ", AMR:" + cls[Freqs_flgs.get("ExAC_AMR")]
                                + ", FIN" + cls[Freqs_flgs.get("ExAC_FIN")] + ", NFE:" + cls[Freqs_flgs.get("ExAC_NFE")]
                                + ", SAS:" + cls[Freqs_flgs.get("ExAC_SAS")] + ", OTH:" + cls[Freqs_flgs.get("ExAC_OTH")];
                        String Freq_gnomAD_exome_POPs = "AFR:" + cls[Freqs_flgs.get("AF_afr")] + ", AMR:" + cls[Freqs_flgs.get("AF_amr")]
                                + ", ASJ" + cls[Freqs_flgs.get("AF_asj")] + ", EAS" + cls[Freqs_flgs.get("AF_eas")]
                                + ", FIN" + cls[Freqs_flgs.get("AF_fin")] + ", NFE:" + cls[Freqs_flgs.get("AF_nfe")]
                                + ", SAS:" + cls[Freqs_flgs.get("AF_sas")] + ", OTH:" + cls[Freqs_flgs.get("AF_oth")];
                        String OMIM = ".";
                        String pheno_MIM;

                        String mim2 = null;
                        String ge = cls[Funcanno_flgs.get("Gene")];
// System.out.println(Funcanno_flgs.get("ExonicFunc.refGene"));
                        String[] gee = ge.split(", ");
                        for (String ggee : gee) {
                            if (mim2gene_dict2.containsKey(ggee)) {
                                mim2 = mim2gene_dict2.get(ggee).toString();
                            } else {
                                mim2 = ".";
                            }
                        }

                        String mim1 = null;
                        String ens = cls[Funcanno_flgs.get("Gene.ensGene")];
                        String[] enss = ens.split(";");
                        for (String eenss : enss) {
                            if (mim2gene_dict.containsKey(eenss)) {
                                mim1 = mim2gene_dict.get(eenss).toString();
                            } else {
                                mim1 = ".";
                            }
                        }

                        if ((mim1 != null) && (!mim1.equals("."))) {
                            OMIM = mim1;
                        }

                        if ((mim2 != null) && (!mim2.equals("."))) {
                            OMIM = mim2;
                        }

                        if (mim_pheno_dict.get(OMIM) == null) {
                            pheno_MIM = ".";
                        } else {
                            pheno_MIM = mim_pheno_dict.get(OMIM).toString();
                        }

                        String orpha_details = "";
                        String orpha_s = " ";
                        String ort3;
                        for (String ort2 : pheno_MIM.split(";")) {
                            if (mim_orpha_dict.get(ort2) == null) {
                                ort3 = ".";
                            } else {
                                ort3 = mim_orpha_dict.get(ort2).toString();
                            }

                            if (!ort3.equals(".")) {
                                orpha_s = ort3 + orpha_s;
                            }
                        }

                        for (String ort4 : orpha_s.split(";")) {
                            if (ort4.length() > 0) {
                                orpha_details = orpha_details + orpha_dict.get(ort4) + "~";
                            }
                        }

                        if (orpha_s.equals("")) {
                            orpha_s = ".";
                        }

                        if (orpha_details.equals("")) {
                            orpha_s = ".";
                        }

                        String thera_out;
                        thera_out = Check_Thera_out(line, Funcanno_flgs);
                        if (thera_out.equals("")) {
                            thera_out = ".";
                        }

                        if (cls.length > Funcanno_flgs.get("Otherinfo")) {
                            writer.write(cls[Allels_flgs.get("Chr")] + "\t" + cls[Allels_flgs.get("Start")] + "\t"
                                    + cls[Allels_flgs.get("End")] + "\t" + cls[Allels_flgs.get("Ref")] + "\t" + cls[Allels_flgs.get("Alt")] + "\t"
                                    + cls[Funcanno_flgs.get("Gene")] + "\t" + cls[Funcanno_flgs.get("Func.refGene")] + "\t"
                                    + cls[Funcanno_flgs.get("ExonicFunc.refGene")] + "\t"
                                    + cls[Funcanno_flgs.get("Gene.ensGene")] + "\t"
                                    + cls[Funcanno_flgs.get("avsnp150")] + "\t clinvar: " + clinvar_bp + " \t VIC: " + vic_bp + " \t" + cls[Funcanno_flgs.get("AAChange.ensGene")] + "\t"
                                    + cls[Funcanno_flgs.get("AAChange.refGene")] + "\t" + cls[Funcanno_flgs.get("AAChange.knownGene")] + "\t"
                                    + cls[Freqs_flgs.get("ExAC_ALL")] + "\t" + cls[Freqs_flgs.get("esp6500siv2_all")] + "\t"
                                    + cls[Freqs_flgs.get("1000g2015aug_all")] + "\t" + cls[Freqs_flgs.get("AF")] + "\t" + cls[Funcanno_flgs.get("FATHMM_pred")] + "\t"
                                    + cls[Funcanno_flgs.get("CADD_phred")] + "\t" + cls[Funcanno_flgs.get("SIFT_score")] + "\t"
                                    + cls[Funcanno_flgs.get("GERP++_RS")] + "\t" + cls[Funcanno_flgs.get("Polyphen2_HDIV_pred")] + "\t"
                                    + cls[Funcanno_flgs.get("dbscSNV_ADA_SCORE")] + "\t" + cls[Funcanno_flgs.get("dbscSNV_RF_SCORE")] + "\t"
                                    + cls[Funcanno_flgs.get("Interpro_domain")] + "\t" + cls[Funcanno_flgs.get("MetaSVM_score")] + "\t" + cls[Funcanno_flgs.get("MutationTaster_pred")] + "\t"
                                    + Freq_ExAC_POPs + "\t" + Freq_gnomAD_exome_POPs + "\t" + OMIM + "\t" + pheno_MIM + "\t" + orpha_s + "\t" + orpha_details + "\t" + thera_out + "\t" + civicdd + "\t" + cls[Funcanno_flgs.get("Otherinfo")] + "\n");
                        } else {
                            writer.write(cls[Allels_flgs.get("Chr")] + "\t" + cls[Allels_flgs.get("Start")] + "\t"
                                    + cls[Allels_flgs.get("End")] + "\t" + cls[Allels_flgs.get("Ref")] + "\t" + cls[Allels_flgs.get("Alt")] + "\t"
                                    + cls[Funcanno_flgs.get("Gene")] + "\t" + cls[Funcanno_flgs.get("Func.refGene")] + "\t"
                                    + cls[Funcanno_flgs.get("ExonicFunc.refGene")] + "\t" + cls[Funcanno_flgs.get("Gene.ensGene")] + "\t"
                                    + cls[Funcanno_flgs.get("avsnp150")] + "\t clinvar: " + clinvar_bp + " \t VIC: " + vic_bp + " \t" + cls[Funcanno_flgs.get("AAChange.ensGene")] + "\t"
                                    + cls[Funcanno_flgs.get("AAChange.refGene")] + "\t" + cls[Funcanno_flgs.get("AAChange.knownGene")] + "\t"
                                    + cls[Freqs_flgs.get("ExAC_ALL")] + "\t" + cls[Freqs_flgs.get("esp6500siv2_all")] + "\t"
                                    + cls[Freqs_flgs.get("1000g2015aug_all")] + "\t" + cls[Freqs_flgs.get("AF")] + "\t" + cls[Funcanno_flgs.get("FATHMM_pred")] + "\t"
                                    + cls[Funcanno_flgs.get("CADD_phred")] + "\t" + cls[Funcanno_flgs.get("SIFT_score")] + "\t"
                                    + cls[Funcanno_flgs.get("GERP++_RS")] + "\t" + cls[Funcanno_flgs.get("Polyphen2_HDIV_pred")] + "\t"
                                    + cls[Funcanno_flgs.get("dbscSNV_ADA_SCORE")] + "\t" + cls[Funcanno_flgs.get("dbscSNV_RF_SCORE")] + "\t"
                                    + cls[Funcanno_flgs.get("Interpro_domain")] + "\t" + cls[Funcanno_flgs.get("MetaSVM_score")] + "\t" + cls[Funcanno_flgs.get("MutationTaster_pred")] + "\t"
                                    + Freq_ExAC_POPs + "\t" + Freq_gnomAD_exome_POPs + "\t" + OMIM + "\t" + pheno_MIM + "\t" + orpha_s + "\t" + orpha_details + "\t" + thera_out + "\t" + civicdd + "\n");
                        }

                    }
                    line_sum += 1;
                }
            }

            writer.flush();
        }

        return line_sum;
    }
}
