#!/usr/bin/env python3
# -*- coding: utf8 -*-

import os
import re
import sys
import shutil
import pathlib
import libs.alinhamento_poli as alinhamento_poli
import libs.alinhamento_poli_enterobacter as alinhamento_poli_enterobacter
import libs.alinhamento_poli_acineto as alinhamento_poli_acineto
import libs.alinhamento_poli_pseudo as alinhamento_poli_pseudo
import libs.alinhamento_outros_pseudo as alinhamento_outros_pseudo
import libs.alinhamento_outros_kleb as alinhamento_outros_kleb
from libs.tools import _bn, _str, sp_runner, count_kraken_words, get_abricate_result, MongoSaver


sys.path[0:0] = ['/opt/pipeline/lib/']


MLST_result = ''
_d = None
a = ''
b = ''
count2 = 0
docker = []
especieAb = []
especieEc = []
especieKp = []
fasta_outros = ''
fasta_polimixina = ''
identificacao = []
lines = []
lista_acineto = ''
lista_enterobacter = ''
lista_kleb = ''
nearest_sts = ''
preidentificacao = []
resultadoANI = ''
THREADS = "16"  # new one will have 32

sys.argv = sys.argv[1:]  # dict(perllib.Array(sys.argv)[1:])
# use alinhamento_poli_truncation;

# Linha de comando do script perl pipeline_melise_output_gal.pl  <caminho do diretorio Output, com resultados de spades/unicycler e prokka ex:/home/melise/Output_run18-06.02.21> <nome da amostra na pasta Output ex:27563_S12> <caminho onde esta instalado o abricate  ex:./abricate/bin/abricate>
# <arquivo da tabela excel formato .xls onde vão ser impressos os resultados> <diretorio onde esta instalado o kmer-db  ex:/home/melise/kmer-db>  <caminho do diretorio onde esta instalado o mlst ex:/home/melise/identificar_clones> <caminho do DB com as seqs de ptn de resistencia mutacional a
# polimixina ex:/home/melise/resis_poli> <caminho do DB com as seqs de ptn de resistencia mutacional a outros antibioticos ex:/home/melise/outrasMut> <caminho kraken ex:kraken> <caminho unicycle ex:unicycler> <caminho R1> <caminho R2>

# Guardar o caminho do diretorio "Output", com o resultado da montagem, anotacao e blast para todas as amostra. ex /home/melise/Output
caminho1 = sys.argv[0]

# Guardar o nome da amostra com o _S*
sample = sys.argv[1]
# para pegar o numero da amostra antes do _S1
sample1 = _str(sample).split('_')
sample2 = sample1[0]
print(f"Sample: {_bn(sample2)}")

# criar um diretório para a amostra
os.makedirs(f"{_bn(caminho1)}/{_bn(sample)}", exist_ok=True)

# caminho para onde esta instalado o abricate ex: ./abricate/bin/abricate

caminho_abricate = sys.argv[2]

# Caminho para o output onde esta a tabela que sera colocado os resultados
caminho_output = sys.argv[3]

# entrar com o caminho da pastar onde esta instalado o kmer-db  /home/melise/kmer-db
kmerdb_install = sys.argv[4]

# entrar com o caminho da pastar onde esta instalado o mlst. ex: /home/melise/identificar_clones
mlst_install = sys.argv[5]

# Caminho para banco de dados com sequências de genes de R a polimixina
db_polimixina = sys.argv[6]

# Caminho para banco de dados com sequências de genes de R mutacionais a outros antibioticos
db_outrosMut = sys.argv[7]

# entrar com o caminho da pastar onde esta instalado o kraken2 ex: kraken2
kraken2_install = sys.argv[8]

# entrar com o caminho do unicycler
unicycler = sys.argv[9]

# criar um diretorio para unicycler
# my @unicycler_dir = ("mkdir","$caminho1/$sample/unicycler");
#		system(@unicycler_dir) == 0
#			or die "system @unicycler_dir failes: $?";

print('Parametros: ')
print(f"caminho: {_bn(caminho1)} ")
print(f"Sample: {_bn(sample)} ")
print(f"SAmple2: {_bn(sample2)} ")
print(f"camino abricate: {_bn(caminho_abricate)} ")
print(f"camino abricate caminho_output: {_bn(caminho_output)} ")
print(f"camino abricate kmerdb_install: {_bn(kmerdb_install)} ")
print(f"mlst install: {_bn(mlst_install)}  ")
print(f"db polimixina: {_bn(db_polimixina)}  ")
print(f"db outros mut: {_bn(db_outrosMut)}  ")
print(f"kraken2_install: {_bn(kraken2_install)}  ")
print(f"unicycler: {_bn(unicycler)} ")

R1 = sys.argv[10]
R2 = sys.argv[11]

mongo_saver = MongoSaver(caminho_output, sample2)  # mongodb instance saver

##################################################
# rodar unicycler
unicycler_exe = " ".join(['unicycler', '-1', f"{_bn(R1)}", '-2', f"{_bn(R2)}", '-o', f"{_bn(caminho1)}/{_bn(sample)}/unicycler", '--min_fasta_length', '500', '--mode', 'conservative', '-t', THREADS, '--spades_path', '/opt/SPAdes-3.13.0-Linux/bin/spades.py'])
sp_runner(unicycler_exe)

# arquivo assembly.fasta da amostra

montagem = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta"

###############################################################################
# rodar prokka
prokka_exe = " ".join(['prokka', '--outdir', f"{_bn(caminho1)}/{_bn(sample)}/prokka", '--prefix', 'genome', f"{montagem}", '--force', "--cpus", "0"])  # --cpus 0 is ALL
sp_runner(prokka_exe)

# EXCLUSIVO DO PIPELINE OUTPUT

# variavel para guardar os tipos de resultados que devem ser indicados
# ex: checkm; especie; contigs; resfinder; VFDB; plasmid; mlst; mutacoes_poli; mutacoes_outra

tipo_de_resultado = None
# o que imprimir desse resultado
imprimir = None

#####################################################################################
# PARA IMPRIMIR O RESULTADO EM UM ARQUIVO TXT PARA GAL

gal_file = open('resultado_gal.txt', mode="a", encoding='utf-8')

##printar no arquivo final o nome da amostra na primeira linha
gal_file.write(f"\nAmostra {_bn(sample)}\nResultados relevantes do sequenciamento do genoma total (WGS):\n")

####################################################################################
# SILENCIADO NO PIPELINE OUTPUT

# my $output = "$caminho_output/$sample.xlsx";
# open (OUT2, ">> $output") or die "Nao foi possivel abrir output\n";

##printar no arquivo final o nome da amostra na primeira coluna
# print OUT2 "$sample\t";

#####################################################################################
# SE NAO QUISER CHECKM

# Rodar o CheckM para saber qualidade da corrida
# Copiar o arquivo assembly.fasta para a pasta do CheckM checkM_bins
os.makedirs("checkM_bins", exist_ok=True)
shutil.copy(os.path.join(".", f"{montagem}"), os.path.join(".", 'checkM_bins'))

# rodar o CheckM
checkM = " ".join(['checkm', 'lineage_wf', '-x', 'fasta', 'checkM_bins', 'checkM_bins', "--threads", THREADS, "--pplacer_threads", THREADS])
sp_runner(checkM)

checkM_qa = " ".join(['checkm', 'qa', '-o', '2', '-f', f"checkM_bins/{_bn(sample)}_resultados", '--tab_table', 'checkM_bins/lineage.ms', 'checkM_bins', "--threads", THREADS])
sp_runner(checkM_qa)

# apagar arquivos gerados, deixando apenas resultados
shutil.rmtree('checkM_bins/bins', ignore_errors=True)
shutil.rmtree('checkM_bins/storage', ignore_errors=True)
pathlib.Path('checkM_bins/assembly.fasta').unlink(missing_ok=True)
pathlib.Path('checkM_bins/lineage.ms').unlink(missing_ok=True)
pathlib.Path('checkM_bins/checkm.log').unlink(missing_ok=True)

# Ler o arquivo do resultado e imprimir o que interessa na tabela
# pegar o resultado da contaminacao

contaminacao = 0

print('Salvando resultado no mongo relatorios')

genome_size = None

with open(f"checkM_bins/{_bn(sample)}_resultados") as IN_check:
    next(IN_check)  # ignore header
    for row in IN_check:
        # remove \n of the line end
        row = row.rstrip("\n")
        # separar as colunas do arquivo em elementos em um array
        lines = row.split("\t")
        # print "$lines[2]\n";
        # printar na tabela OUTPUT as informacoes de qualidade que interessam EXCLUSIVO TABELA OUTPUT
        # print OUT2 "$lines[5]\t$lines[6]\t$lines[8]\t";
        genome_size = lines[8]
        mongo_saver.save('checkm_1', lines[5])
        mongo_saver.save('checkm_2', lines[6])
        mongo_saver.save('checkm_3', lines[8])
        mongo_saver.save('checkm_4', lines[11])  # contigs
        contaminacao = lines[6]

mongo_saver.save('sample', sample)
#########################################################################################################################
# Identificar especie usando o kraken

print('rodar o kraken')
kraken = " ".join([f"{_bn(kraken2_install)}/kraken2", '--db', f"{_bn(kraken2_install)}/minikraken2_v2_8GB_201904_UPDATE", '--use-names', '--paired', f"{_bn(R1)}", f"{_bn(R2)}", '--output', 'out_kraken', "--threads", THREADS])
sp_runner(kraken)

print("splitting output into %s equal files" % THREADS)
preffix = "krk"
splitter = " ".join(["split", "--numeric-suffixes=1", "-n", f"l/{THREADS}", "out_kraken", preffix])
sp_runner(splitter)
ordenado = count_kraken_words(int(THREADS), preffix)

# contar qts vezes aparece cada especie
repeticoes = ordenado.most_common(2)
maior_repeticao = repeticoes[0][0]
segunda_repeticao = repeticoes[1][0]

# print "$maior_repeticao\n$segunda_repeticao\n";

# onde sera guardada o nome da especie
check_especies = maior_repeticao

# apagar o arquivo do resultado
# my @rm2 = ("rm", "out_kraken");

# system(@rm2) == 0
#        or die "system @rm2 failes: $?";

# colocar só o genero e a especie, descartando qualquer outra informação dessa coluna do kraken
identificar_especie = ''  # mod 11.05.22
genero = ''  # mod 11.05.22
especie = ''  # mod 11.05.22
# juntar genero e especie para mlst
# resultado_final_especie = ''  # mod 11.05.22
# resultado que sera impresso
printar_especies = ''  # mod 11.05.22
# o que usar para mlst
especie_mlst = ''  # mod 11.05.22

if (re.findall(re.compile(r'\w+\s\w+', re.I), check_especies)):  # mod 11.05.22
    check_especies = check_especies.strip()
    genero, especie = check_especies.split(" ")
    # print "$genero\n$especie\n";
    # Associar o nome da especie ao banco do mlst e gravar o nome que sera dado como resultado final
    # resultado_final_especie = f"{genero}{especie}"
    # $printar_especies = $resultado_final_especie;
    # $especie_mlst = "";
else:
    printar_especies = check_especies
    genero = check_especies
# mod ate aqui 20.05.22

#######################################################################################################################
# Sequencia para verificar mutacoes pontuais usando subrotinas proprias

# guardar o resultado das mutacoes para polimixina

result2 = []
# guardar o resultado das mutacoes para outros antibioticos
result3 = []
# guardar resultados dos fatores de virulencia
vfdb = []

resultado_final_especie = f"{genero}{especie}".lower()
print(f"resultado_final_especie: {resultado_final_especie}")
if resultado_final_especie == 'pseudomonasaeruginosa':
    especie_mlst = 'paeruginosa'
    printar_especies = 'Pseudomonas aeruginosa'
    fasta_polimixina = f"{_bn(db_polimixina)}/proteins_pseudo_poli.fasta"
    result2 = alinhamento_poli_pseudo.poli_pseudo(montagem, fasta_polimixina, sample, THREADS)
    fasta_outros = f"{_bn(db_outrosMut)}/proteins_outrasMut_pseudo.fasta"
    result3 = alinhamento_outros_pseudo.outros_pseudo(montagem, fasta_outros, sample, THREADS)
elif resultado_final_especie == 'escherichiacoli':
    especie_mlst = 'ecoli'
    printar_especies = 'Escherichia coli'
elif resultado_final_especie == 'staphylococcusaureus':
    especie_mlst = 'saureus'
    printar_especies = 'Staphylococcus aureus'
elif resultado_final_especie == 'pseudomonasputida':
    especie_mlst = 'pputida'
    printar_especies = 'pseudomonas putida'
elif resultado_final_especie == 'Listeriamonocytogenes':  # modificado 19.11.21
    especie_mlst = 'lmonocytogenes'
    printar_especies = 'Listeria monocytogenes'
elif resultado_final_especie == 'enterococcusfaecalis':
    especie_mlst = 'efaecalis'
    printar_especies = 'Enterococcus faecalis'
elif resultado_final_especie == 'klebsiellaoxytoca':
    especie_mlst = 'koxytoca'
    printar_especies = 'Klebsiella oxytoca'
elif resultado_final_especie == 'enterococcusfaecium':
    especie_mlst = 'efaecium'
    printar_especies = 'Enterococcus faecium'
elif resultado_final_especie == 'serratiamarcescens':
    especie_mlst = 'Nao disponivel'
    printar_especies = 'Serratia marcescens'
elif resultado_final_especie == 'providenciastuartii':
    especie_mlst = 'Nao disponivel'
    printar_especies = 'Providencia stuartii'
elif resultado_final_especie in ('klebsiellapneumoniae', 'acinetobacterbaumannii',
                                 "enterobactercloacae", "enterobacterhormaechei", "enterobacterasburiae",
                                 "enterobacterkobei", "enterobacterroggenkampii", "enterobacterludwigii"):
    lista = ""
    fastANI_txt = "Para FastANI"
    if resultado_final_especie == 'klebsiellapneumoniae':
        especie_mlst = 'kpneumoniae'
        # $printar_especies = "Klebsiella pneumoniae";
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_kleb_poli.fasta"
        result2 = alinhamento_poli.poli(montagem, fasta_polimixina, sample, THREADS)
        # @result2 = &alinhamento_poli_truncation::poli($montagem,$fasta_polimixina,$sample);
        fasta_outros = f"{_bn(db_outrosMut)}/proteins_outrasMut_kleb.fasta"
        result3 = alinhamento_outros_kleb.outros_kleb(montagem, fasta_outros, sample, THREADS)
        lista = '/opt/genomas_enterobacter/kleb_database/lista-kleb'  # CAMBIAR
    elif resultado_final_especie == "acinetobacterbaumannii":
        especie_mlst = 'abaumannii_2'
        # $printar_especies = "Acinetobacter baumannii";
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_acineto_poli.fasta"
        result2 = alinhamento_poli_acineto.poli_acineto(montagem, fasta_polimixina, sample, THREADS)
        lista = '/opt/genomas_enterobacter/fastANI_acineto/list-acineto'  # CAMBIAR
    elif resultado_final_especie in ("enterobactercloacae", "enterobacterhormaechei", "enterobacterasburiae",
                                     "enterobacterkobei", "enterobacterroggenkampii", "enterobacterludwigii"):
        especie_mlst = "ecloacae"
        lista = '/opt/genomas_enterobacter/fastANI/list_entero'  # CAMBIAR
        fastANI_txt = 'Rodar fastANI para subespecie'

    print(fastANI_txt)
    # Abrir o arquivo lista

    fastani = " ".join(['/opt/FastANI/fastANI', '-q', f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", '--rl', f"{lista}", '-o', f"{_bn(sample)}_out-fastANI", "--threads", THREADS])
    sp_runner(fastani)

    # abrir output
    # Abrir o arquivo do output de distancia
    # array para guardar especies

    print('resultado do fastANI')
    with open(f"{_bn(sample)}_out-fastANI", "r") as IN7:
        especiE = IN7.readline().rstrip("\n").split('\t')  # first line only
        preidentificacao = especiE[1].split("/")
        identificacao = preidentificacao[-1].split(".")
        printar_especies = identificacao[0]

    if re.search(r'Enterobacter_cloacae_subsp_cloacae', printar_especies, re.IGNORECASE):
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_Ecloacae_poli.fasta"
        result2 = alinhamento_poli_enterobacter.poli_enterobacter(montagem, fasta_polimixina, sample, THREADS)
else:
    # mod 20.05.22
    printar_especies = f"{genero} {especie}"  # mod 10.05.22
    especie_mlst = 'Nao disponivel'  # mod 26.08.22
# mod 20.05.22

# print "$especie_mlst ";

print('contaminacao...')
# printar no arquivo final o nome da especie
if float(contaminacao) <= 10.:
    # print OUT2 "$printar_especies\t";
    mongo_saver.save('especie', printar_especies)
    # para o gal
    gal_file.write(f"Espécie identificada: {_bn(printar_especies)}\n")
else:
    # print OUT2 "$printar_especies\t";
    imprimir = f"{maior_repeticao} {_bn(repeticoes[0][1])} {segunda_repeticao} {_bn(repeticoes[1][1])}"
    mongo_saver.save('especie', imprimir)
    # para o gal
    gal_file.write(f"Espécie: CONTAMINAÇÃO {_bn(imprimir)}\n")

# else {
#        	print OUT2 "$maior_repeticao $count_ordenado2{$maior_repeticao} $segunda_repeticao $count_ordenado2{$segunda_repeticao}\t";
# }

###############################################################################################
# Rodar ABRICATE
# Para resistencia usando o ResFinder (porque so tem resistencia adquirida)
abricante_out = f"{_bn(sample)}_outAbricateRes"
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "resfinder", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", abricante_out, '--threads', THREADS])
sp_runner(abricante_exe)
selected = get_abricate_result(abricante_out, resistencia=True)
select_imprimir = []

# criar um @ para cada classe de antibioticos
genes = []

# print gal
gal_file.write('Genes de resistência encontrados no WGS:\n')

for n_l in selected:
    # print "$n\n";
    # separar as colunas do arquivo em elementos de um array
    lines_blast = n_l.split("\t")
    # concatenar os resultado
    out_blast = f"{lines_blast[5]} (ID:{lines_blast[10]} COV_Q:{lines_blast[9]} COV_DB:{lines_blast[6]})"
    select_imprimir.append(out_blast)
    # imprimir no arquivo do gal
    if re.match(r'.*(blaKPC|blaNDM|blaVIM|blaIMP|blaSPM|blaOXA-23|blaOXA-24|blaOXA-25|blaOXA-26|blaOXA-27|blaOXA-48|blaOXA-51|blaOXA-58|blaOXA-64|blaOXA-65|blaOXA-69|blaOXA-90|blaOXA-72|blaOXA-98|blaOXA-116|blaOXA-117|blaOXA-160|blaOXA-175|blaOXA-176|blaOXA-253|blaOXA-343).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (carbapenemase)\n";
        genes.append(lines_blast[5] + " (carbapenemase)")
    elif re.match(r'.*(blaTEM|blaSHV|blaADC|blaCTX-M|blaGES|blaOXA-(?!23|24|25|26|27|48|51|58|64|65|69|72|98|90|116|117|160|175|176|253|343)).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (ESBL)\n";
        genes.append(lines_blast[5] + " (ESBL)")
    elif re.match(r"(aac\(6\'\)-Ib-cr).*", lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a aminoglicosídeos e fluoroquinolonas)\n";
        genes.append(f"{lines_blast[5]} (resistance to aminoglycosides and fluoroquinolones)")
    elif re.match(r'(aph|aac|rmt|aad).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a aminoglicosídeos)\n";
        genes.append(f"{lines_blast[5]} (resistance to aminoglycosides)")
    elif re.match(r'(cat|cml|cmx|floR).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia ao cloranfenicol)\n";
        genes.append(f"{lines_blast[5]} (resistance to chloramphenicol)")
    elif re.match(r'(qnr|oqx).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a fluoroquinolonas)\n";
        genes.append(f"{lines_blast[5]} (resistance to fluoroquinolones)")
    elif re.match(r'sul.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a sulfonamidas)\n";
        genes.append(f"{lines_blast[5]} (resistance to sulfonamidas)")
    elif re.match(r'dfrA.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a trimetoprim)\n";
        genes.append(f"{lines_blast[5]} (resistance to trimetoprim)")
    elif re.match(r'tet.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a tetraciclina)\n";
        genes.append(f"{lines_blast[5]} (resistance to tetracycline)")
    elif re.match(r'ere.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a eritromicina)\n";
        genes.append(f"{lines_blast[5]} (resistance to eritromicina)")
    elif re.match(r'erm.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a lincosamidas, macrolideos e estreptograminas)\n";
        genes.append(f"{lines_blast[5]} (resistance to lincosamides, macrolides and streptogramins)")
    elif re.match(r'ARR.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a rifampicina)\n";
        genes.append(f"{lines_blast[5]} (resistance to rifampicin)")
    elif re.match(r'(mph|msr).*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a macrolideos)\n";
        genes.append(f"{lines_blast[5]} (resistance to macrolides)")
    elif re.match(r'.*Van.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
        genes.append(f"{lines_blast[5]} (resistance to vancomycin)")
    elif re.match(r'.*lsa.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
        genes.append(f"{lines_blast[5]} (resistance to clindamycin)")
    elif re.match(r'.*mcr.*', lines_blast[5], re.I):
        # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
        genes.append(f"{lines_blast[5]} (resistance to polymyxin)")  # COLOCAR O NOME EM INGLES
    else:
        # print OUT2 "$lines_blast[5]\n";
        genes.append(f"{lines_blast[5]}")

# imprimir resultados com a classe do antimicrobiano

mongo_saver.save("gene", "<br>".join(genes))
mongo_saver.save("resfinder", "<br>".join(select_imprimir))

################################################################################################
# Rodar abricate para VFDB (Virulence factor)
abricante_out = f"{_bn(sample)}_outAbricateVFDB"
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "vfdb", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", abricante_out, '--threads', THREADS])
sp_runner(abricante_exe)
selected = get_abricate_result(abricante_out)
select_imprimir = []

# ler o array selected
for n_l in selected:
    # print "$n\n";
    # separar as colunas do arquivo em elementos de um array
    lines_blast = n_l.split("\t")
    # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
    out_blast = f"{lines_blast[1]}: {lines_blast[5]} {lines_blast[13]} ID:{lines_blast[10]} COV_Q:{lines_blast[9]} COV_DB:{lines_blast[6]}| "
    select_imprimir.append(out_blast)

mongo_saver.save('VFDB', "<br>".join(select_imprimir))
pathlib.Path(abricante_out).unlink(missing_ok=True)

#########################################################################################################
# Rodar abricate para PlasmidFinder
abricante_out = f"{_bn(sample)}_outAbricatePlasmid"
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "plasmidfinder", f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", ">", abricante_out, '--threads', THREADS])
sp_runner(abricante_exe)
selected = get_abricate_result(abricante_out)
select_imprimir = []

# ler o array selected
imprimir = 'Not found'
if not selected:
    # print OUT2 "Nao encontrado\t";
    mongo_saver.save('plasmid', imprimir)
    # para o gal
else:
    imprimir = ""
    for n_l in selected:
        # print "$n\n";
        # separar as colunas do arquivo em elementos de um array
        lines_blast = n_l.split("\t")
        # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
        out_blast = lines_blast[5] + 'ID' + ':' + lines_blast[10] + ' ' + 'COV_Q:' + lines_blast[9] + ' ' + 'COV_DB:' + lines_blast[6] + '|' + ' '
        select_imprimir.append(out_blast)
        imprimir += f"\n{lines_blast[5]}"
gal_file.write(f"Plasmídeos encontrados:{_bn(imprimir)}\n")

mongo_saver.save("plasmid", "<br>".join(select_imprimir))
pathlib.Path(abricante_out).unlink(missing_ok=True)

################################################################################################

print(f"Rodar o MLST {especie_mlst}")

MLST_result = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/data.json"
# se nao tem mlst disponivel, ai tem que avisar
if (especie_mlst == 'Nao disponivel') or (especie_mlst == ''):  # mod 26-08-22
    # print OUT2 "Nao disponivel\t";
    imprimir = 'Not available for this species'  # mod 26.08.22
    mongo_saver.save('mlst', imprimir)
    # para o gal
    gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)\n")
else:
    # mod 26-08-22
    docker = " ".join(['docker', 'run', '--rm', '-i', '-v', f"{_bn(mlst_install)}/mlst_db:/database", '-v', f"{_bn(caminho1)}/{_bn(sample)}/unicycler:/workdir", 'mlst', '-i', 'assembly.fasta', '-o', '.', '-s', f"{especie_mlst}"])
    # rodar o mlst
    sp_runner(docker)
# mod

ST = None
print('ler o resultado do mlst')
mlst_json = pathlib.Path(MLST_result)
mlst_json.touch(exist_ok=True)  # will create file, if it exists will do nothing

with open(mlst_json, "r") as IN3:
    line = IN3.readline().rstrip("\n")  # single line file
    a = re.search(r'.*sequence_type":\s"(\d{1,4})".*', line, re.IGNORECASE)
    b = re.search(r'.*sequence_type":\s"(\d*!,\d*!)".*', line, re.IGNORECASE)
    c = re.search(r'.*sequence_type":\s"(\d{1,4}\*)".*', line, re.IGNORECASE)
    for m in (a, b, c):
        if not m:
            continue
        ST = m.group(1)
        # print OUT2 "$ST\t";
        imprimir = ST
        mongo_saver.save('mlst', imprimir)
        # para o gal
        gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)\n")
    m = re.search(r'nearest_sts":\s"((\d*,)*\d*)".*', line, re.IGNORECASE)
    if m:
        nearest_sts = m.group(1)
        if nearest_sts:
            # print OUT2 "Nearest $nearest_sts\t";
            imprimir = f"Nearest {nearest_sts}"
            mongo_saver.save('mlst', imprimir)
            # para o gal
            gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)\n")
    m = re.search(r'.*sequence_type":\s"(Unknown)".*', line, re.IGNORECASE)
    if m:
        ST = m.group(1)
        # print OUT2 "Unknown\t";
        imprimir = 'Unknown'
        mongo_saver.save('mlst', imprimir)
        # para o gal
        gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)\n")

mongo_saver.save('mutacoes_poli', "<br>".join(result2))
gal_file.write("Mutações polimixina: %s" % "<br>".join(result2))
mongo_saver.save('mutacoes_outras', "<br>".join(result3))

######################################################################
print('rodar coverage')

# figuring out if file is compressed or not
catcmd = "cat"
res = sp_runner(f"file {_bn(R1)}", pipeout=True)
if res and str(res).find("gzip compressed") > -1:
    catcmd = "zcat"

zcat = " ".join([f"echo $({catcmd} {_bn(R1)} | wc -l)/4 | bc"])
res_r1 = sp_runner(zcat, pipeout=True)
n_reads1 = res_r1.decode("utf-8").rstrip("\n")

# o mesmo para o arquivo R2
zcat2 = " ".join([f"echo $({catcmd} {_bn(R2)} | wc -l)/4 | bc"])
res_r2 = sp_runner(zcat2, pipeout=True)
n_reads2 = res_r2.decode("utf-8").rstrip("\n")

soma_reads = (float(n_reads1) + float(n_reads2))

###calcular tamanho medio das reads, vou usar só as R1 como base
zcat3 = " ".join([f"{catcmd} {_bn(R1)} | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}'"])
res_avg = sp_runner(zcat3, pipeout=True)
average_length2 = res_avg.decode("utf-8").rstrip("\n")

gal_file.close()

coverage = (float(average_length2) * soma_reads) / float(genome_size)

mongo_saver.save('coverage', coverage)
