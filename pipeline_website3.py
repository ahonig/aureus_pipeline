#!/usr/bin/env python3
# -*- coding: utf8 -*-

import os
import re
import sys
import shutil
import pathlib
import subprocess
import libs.alinhamento_poli as alinhamento_poli
import libs.alinhamento_poli_enterobacter as alinhamento_poli_enterobacter
import libs.alinhamento_poli_acineto as alinhamento_poli_acineto
import libs.alinhamento_poli_pseudo as alinhamento_poli_pseudo
import libs.alinhamento_outros_pseudo as alinhamento_outros_pseudo
import libs.alinhamento_outros_kleb as alinhamento_outros_kleb
import libs.save_result_mongo as save_result_mongo

_bn = lambda s: '' if s is None else s
_str = lambda s: '' if s is None else str(s)
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
print(f"caminho1: {_bn(caminho1)} ")
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

##################################################
# rodar unicycler
unicycler_exe = " ".join(['unicycler', '-1', f"{_bn(R1)}", '-2', f"{_bn(R2)}", '-o', f"{_bn(caminho1)}/{_bn(sample)}/unicycler", '--min_fasta_length', '500', '--mode', 'conservative', '-t', THREADS, '--spades_path', '/opt/SPAdes-3.13.0-Linux/bin/spades.py'])
p2 = subprocess.Popen(unicycler_exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_uni = p2.communicate()
if p2.returncode:
    raise Exception(f"system {unicycler_exe} failes: {res_uni[1]}")

# arquivo assembly.fasta da amostra

montagem = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta"

###############################################################################
# rodar prokka
prokka_exe = " ".join(['prokka', '--outdir', f"{_bn(caminho1)}/{_bn(sample)}/prokka", '--prefix', 'genome', f"{montagem}", '--force', "--cpus", "0"])  # --cpus 0 is ALL
p2 = subprocess.Popen(prokka_exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_prk = p2.communicate()
if p2.returncode:
    raise Exception(f"system {prokka_exe} failes: {res_prk[1]}")

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
gal_file.write(f"\nAmostra {_bn(sample)}\nResultados relevantes do sequenciamento do genoma total (WGS):")

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
shutil.copyfile(os.path.join(".", f"{montagem}"), os.path.join(".", 'checkM_bins'))

# rodar o CheckM

checkM = " ".join(['checkm', 'lineage_wf', '-x', 'fasta', 'checkM_bins', 'checkM_bins', "--threads", THREADS, "--pplacer_threads", THREADS])
p2 = subprocess.Popen(checkM, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_m = p2.communicate()
if p2.returncode:
    raise Exception(f"system {checkM} failes: {res_m[1]}")

checkM_qa = " ".join(['checkm', 'qa', '-o', '2', '-f', f"checkM_bins/{_bn(sample)}_resultados", '--tab_table', 'checkM_bins/lineage.ms', 'checkM_bins', "--threads", THREADS])
p2 = subprocess.Popen(checkM_qa, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_qa = p2.communicate()
if p2.returncode:
    raise Exception(f"system {checkM_qa} failes: {res_qa[1]}")

# apagar arquivos gerados, deixando apenas resultados
shutil.rmtree('checkM_bins/bins')
shutil.rmtree('checkM_bins/storage')
shutil.rmtree('checkM_bins/assembly.fasta')
shutil.rmtree('checkM_bins/lineage.ms')
shutil.rmtree('checkM_bins/checkm.log')

# Ler o arquivo do resultado e imprimir o que interessa na tabela
# pegar o resultado da contaminacao

contaminacao = 0

print('Salvando resultado no mongo relatorios')

genome_size = None

cabecera = 'S'

with open(f"checkM_bins/{_bn(sample)}_resultados") as IN_check:
    for row in IN_check:
        # remove \n of the line end
        row = row.rstrip("\n")
        # separar as colunas do arquivo em elementos em um array
        if cabecera == 'N':
            lines = row.split("\t")
            # print "$lines[2]\n";
            # printar na tabela OUTPUT as informacoes de qualidade que interessam EXCLUSIVO TABELA OUTPUT
            # print OUT2 "$lines[5]\t$lines[6]\t$lines[8]\t";
            tipo_de_resultado = 'checkm_1'
            imprimir = lines[5]
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            tipo_de_resultado = 'checkm_2'
            imprimir = lines[6]
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            tipo_de_resultado = 'checkm_3'
            imprimir = lines[8]
            genome_size = imprimir
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            tipo_de_resultado = 'checkm_4'  # contigs
            imprimir = lines[11]
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            contaminacao = lines[6]
        cabecera = 'N'

tipo_de_resultado = 'sample'
imprimir = sample
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

#########################################################################################################################
# Identificar especie usando o kraken

# onde sera guardada o nome da especie
check_especies = ''

print('rodar o kraken')
kraken = " ".join([f"{_bn(kraken2_install)}/kraken2", '--db', f"{_bn(kraken2_install)}/minikraken2_v2_8GB_201904_UPDATE", '--use-names', '--paired', f"{_bn(R1)}", f"{_bn(R2)}", '--output', 'out_kraken', "--threads", THREADS])
p2 = subprocess.Popen(kraken, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_krk = p2.communicate()
if p2.returncode:
    raise Exception(f"system {kraken} failes: {res_krk[1]}")

ordenado = []
ordenado2 = []
with open("out_kraken", "r") as IN6:
    for row in IN6:
        row = row.rstrip('\n')
        out_kraken2 = row.split(r'	')
        ordenado.append(out_kraken2[2])

    for n in ordenado:
        a = re.search(r'(\w*\s\w*\s).*\(taxid.*\)', n, re.IGNORECASE)
        b = re.search(r'(\w+\ssp\.).*\(taxid.*\)', n, re.IGNORECASE)
        c = re.search(r'^(\w+)\s\(taxid.*\)', n, re.IGNORECASE)
        for x in (a, b, c):
            if not x:
                continue
            ordenado2.append(x.group(1))

# contar qts vezes aparece cada especie

count_ordenado2 = {}
for item in ordenado2:
    if item not in count_ordenado2:
        count_ordenado2[item] = 0
    count_ordenado2[item] += 1

repeticoes = sorted(count_ordenado2.keys(), key=lambda x: count_ordenado2[x])
maior_repeticao = repeticoes[-1]
segunda_repeticao = ""  # repeticoes[-2]

# print "$maior_repeticao\n$segunda_repeticao\n";

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
        lista = '/opt/genomas_enterobacter/fastANI/list_entero'  # CAMBIAR
        fastANI_txt = 'Rodar fastANI para subespecie'

    print(fastANI_txt)
    # Abrir o arquivo lista

    fastani = " ".join(['/opt/FastANI/fastANI', '-q', f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", '--rl', f"{lista}", '-o', f"{_bn(sample)}_out-fastANI", "--threads", THREADS])
    p2 = subprocess.Popen(fastani, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    res_fani = p2.communicate()
    if p2.returncode:
        raise Exception(f"system {fastani} failes: {res_fani[1]}")

    # abrir output
    # Abrir o arquivo do output de distancia
    # array para guardar especies
    count2 = 0
    print('resultado do fastANI')
    with open(f"{_bn(sample)}_out-fastANI", "r") as IN7:
        especiE = IN7.readline().rstrip("\n").split('\t')  # first line only
        preidentificacao = especiE[1].split("/")
        identificacao = preidentificacao[-1].split("\.")  # TODO: confusion here about the index
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
    tipo_de_resultado = 'especie'
    imprimir = printar_especies
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    gal_file.write(f"Espécie identificada: {_bn(imprimir)}")

if float(contaminacao) > 10.:
    # print OUT2 "$printar_especies\t";
    tipo_de_resultado = 'especie'
    imprimir = f"{maior_repeticao} {_bn(count_ordenado2.get(maior_repeticao, ''))} {segunda_repeticao} {_bn(count_ordenado2.get(segunda_repeticao, ''))}"
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    gal_file.write(f"Espécie: CONTAMINAÇÃO {_bn(imprimir)}")

# else {
#        	print OUT2 "$maior_repeticao $count_ordenado2{$maior_repeticao} $segunda_repeticao $count_ordenado2{$segunda_repeticao}\t";
# }

###############################################################################################
# Rodar ABRICATE
# Para resistencia usando o ResFinder (porque so tem resistencia adquirida)
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "resfinder", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", f"{_bn(sample)}_outAbricateRes", '--threads', THREADS])
p2 = subprocess.Popen(abricante_exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_abri = p2.communicate()
if p2.returncode:
    raise Exception(f"system {abricante_exe} failes: {res_abri[1]}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# criar um array para guardar so os que tiverem identidade e cobertura alta

selected = []
rowQty = 0
with open(f"{_bn(sample)}_outAbricateRes", "r") as IN8:
    for row in IN8:
        rowQty += 1
        if rowQty <= 4:
            continue
        # remove \n of the line end
        row = row.rstrip("\n")
        # print "$row\n";
        # separar as colunas do arquivo em elementos em um array
        lines = row.split("\t")
        # print "$lines[10]\n";
        # printar no array selected so quem tem identidade maior que 90
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected.append(f"{_bn(row)}")
            print(f"{_bn(row)}")
        # if para o operon de Van no Resfinder

        if re.match(r'Van.*', lines[5], re.I):
            selected.append(f"{_bn(row)}")
            print(f"{_bn(row)}")

# Analisar o array dos contigs
# ler o array selected


select_imprimir = []
out_blast = ''

# criar um @ para cada classe de antibioticos
genes = []

# print gal
gal_file.write('Genes de resistência encontrados no WGS: ')

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

# zerar a variarvel para concatenar
imprimir = ''

for _d in genes:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'gene'  # MOD
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

# zerar a variarvel para concatenar
imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in select_imprimir:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'resfinder'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

################################################################################################
# Rodar abricate para VFDB (Virulence factor)
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "vfdb", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", f"{_bn(sample)}_outAbricateVFDB", '--threads', THREADS])
p2 = subprocess.Popen(abricante_exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_abri2 = p2.communicate()
if p2.returncode:
    raise Exception(f"system {abricante_exe} failes: {res_abri2[1]}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# abrir o resultado do abricate_VFDB
# criar um array para guardar so os que tiverem identidade e cobertura alta

selected3 = []
select_imprimir3 = []
out_blast3 = ''
qty = 0
with open(f"{_bn(sample)}_outAbricateVFDB", "r") as IN10:
    for row in IN10:
        qty += 1
        if qty <= 4:
            continue
        # remove \n of the line end
        row = row.rstrip("\n")
        # print "$row\n";
        # separar as colunas do arquivo em elementos em um array
        lines = row.split("\t")
        # print "$lines[10]\n";
        # printar no array selected so quem tem coverage maior que 90
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected3.append(f"{_bn(row)}")

# ler o array selected
for n_l in selected3:
    # print "$n\n";
    # separar as colunas do arquivo em elementos de um array
    lines_blast = n_l.split("\t")
    # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
    out_blast3 = f"{lines_blast[1]}: {lines_blast[5]} {lines_blast[13]} ID:{lines_blast[10]} COV_Q:{lines_blast[9]} COV_DB:{lines_blast[6]}| "
    select_imprimir3.append(out_blast3)

# zerar a variarvel para concatenar
imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in select_imprimir3:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'VFDB'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

shutil.rmtree(f"{_bn(sample)}_outAbricateVFDB")

#########################################################################################################
# Rodar abricate para PlasmidFinder
abricante_exe = " ".join([f"{_bn(caminho_abricate)}", "--db", "plasmidfinder", f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", ">", f"{_bn(sample)}_outAbricatePlasmid", '--threads', THREADS])
p2 = subprocess.Popen(abricante_exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_abri3 = p2.communicate()
if p2.returncode:
    raise Exception(f"system {abricante_exe} failes: {res_abri3[1]}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# abrir o resultado do abricate_Plasmid
# criar um array para guardar so os que tiverem identidade e cobertura alta
selected4 = []
select_imprimir4 = []
out_blast4 = None
qxty = 0

with open(f"{_bn(sample)}_outAbricatePlasmid", "r") as IN11:
    for row in IN11:
        qxty += 1
        if qxty <= 4:
            continue
        # remove \n of the line end
        row = row.rstrip("\n")
        # print "$row\n";
        # separar as colunas do arquivo em elementos em um array
        lines = row.split("\t")
        # print "$lines[10]\n";
        # printar no array selected so quem tem coverage maior que 90
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected4.append(f"{_bn(row)}")
# ler o array selected

if not selected4:
    # print OUT2 "Nao encontrado\t";
    tipo_de_resultado = 'plasmid'
    imprimir = 'Not found'
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    gal_file.write(f"Plasmídeos encontrados: {_bn(imprimir)}")
else:
    gal_file.write('Plasmídeos encontrados:')
    for n_l in selected4:
        # print "$n\n";
        # separar as colunas do arquivo em elementos de um array
        lines_blast = n_l.split("\t")
        # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
        out_blast4 = lines_blast[5] + 'ID' + ':' + lines_blast[10] + ' ' + 'COV_Q:' + lines_blast[9] + ' ' + 'COV_DB:' + lines_blast[6] + '|' + ' '
        select_imprimir4.append(out_blast4)
        gal_file.write(f"{_bn(lines_blast[5])}")
    # print OUT2 "\t";

# zerar a variarvel para concatenar
imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir se não tiver vazia
if not select_imprimir4:
    # nao fazer nada se estiver vazio
    pass
else:
    for _d in select_imprimir4:
        imprimir = _str(imprimir) + _str(_d)
        imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'plasmid'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

shutil.rmtree(f"{_bn(sample)}_outAbricatePlasmid")

################################################################################################

print(f"Rodar o MLST {especie_mlst}")

MLST_result = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/data.json"
# se nao tem mlst disponivel, ai tem que avisar
if (especie_mlst == 'Nao disponivel') or (especie_mlst == ''):  # mod 26-08-22
    # print OUT2 "Nao disponivel\t";
    tipo_de_resultado = 'mlst'
    imprimir = 'Not available for this species'  # mod 26.08.22
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)")
else:
    # mod 26-08-22
    docker = " ".join(['docker', 'run', '--rm', '-i', '-v', f"{_bn(mlst_install)}/mlst_db:/database", '-v', f"{_bn(caminho1)}/{_bn(sample)}/unicycler:/workdir", 'mlst', '-i', 'assembly.fasta', '-o', '.', '-s', f"{especie_mlst}"])

    # rodar o mlst
    p2 = subprocess.Popen(docker, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    res_docker = p2.communicate()
    if p2.returncode:
        raise Exception(f"system {docker} failes: {res_docker[1]}")
# mod

ST = None
print('ler o resultado do mlst')
mlst_file = pathlib.Path(MLST_result)
mlst_file.touch(exist_ok=True)  # will create file, if it exists will do nothing

with open(mlst_file, "r") as IN3:
    for row3 in IN3:
        row3 = row3.rstrip("\n")
        a = re.search(r'.*sequence_type":\s"(\d{1,4})".*', row3, re.IGNORECASE)
        b = re.search(r'.*sequence_type":\s"(\d*!,\d*!)".*', row3, re.IGNORECASE)
        c = re.search(r'.*sequence_type":\s"(\d{1,4}\*)".*', row3, re.IGNORECASE)
        for m in (a, b, c):
            if not m:
                continue
            ST = m.group(1)
            # print OUT2 "$ST\t";
            tipo_de_resultado = 'mlst'
            imprimir = ST
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            # para o gal
            gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)")
        m = re.match(r'nearest_sts":\s"((\d*,)*\d*)".*', row3, re.IGNORECASE)
        if m:
            nearest_sts = m.group(1)
            # print OUT2 "Nearest $nearest_sts\t";
            tipo_de_resultado = 'mlst'
            imprimir = f"Nearest {nearest_sts}"
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            # para o gal
            gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)")
        m = re.match(r'.*sequence_type":\s"(Unknown)".*', row3, re.IGNORECASE)
        if m:
            ST = m.group(1)
            # print OUT2 "Unknown\t";
            tipo_de_resultado = 'mlst'
            imprimir = 'Unknown'
            save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
            # para o gal
            gal_file.write(f"Clone ST {_bn(imprimir)} (determinado por MLST)")

# zerar a variarvel para concatenar

imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in result2:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'mutacoes_poli'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
# para o gal
gal_file.write(f"Mutações polimixina: {_bn(imprimir)}\n")

# zerar a variarvel para concatenar
imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in result3:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

# print OUT2 "@result3\t";

tipo_de_resultado = 'mutacoes_outras'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

######################################################################
print('rodar coverage')

# figuring out if file is compressed or not
catcmd = "cat"
ftype = f"file {_bn(R1)}"
p2 = subprocess.Popen(ftype, stdout=subprocess.PIPE, shell=True)
res_neko = p2.communicate()
if res_neko and str(res_neko[0]).find("gzip compressed") > -1:
    catcmd = "zcat"

zcat = " ".join([f"echo $({catcmd} {_bn(R1)} | wc -l)/4 | bc > out_R1 "])
p2 = subprocess.Popen(zcat, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_r1 = p2.communicate()
if p2.returncode:
    raise Exception(f"system {zcat} failes: {res_r1[1]}")

n_reads1 = .0
# guardar o numero de reads
with open('out_R1', "r") as IN:
    for row in IN:
        row = row.rstrip("\n")
        n_reads1 = row

# o mesmo para o arquivo R2
zcat2 = " ".join([f"echo $({catcmd} {_bn(R2)} | wc -l)/4 | bc > out_R2"])
p2 = subprocess.Popen(zcat2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_r2 = p2.communicate()
if p2.returncode:
    raise Exception(f"system {zcat2} failes: {res_r2[1]}")

n_reads2 = .0
# guardar o numero de reads
with open('out_R2', "r") as IN:
    for row in IN:
        row = row.rstrip("\n")
        n_reads2 = row

soma_reads = (float(n_reads1) + float(n_reads2))

###calcular tamanho medio das reads, vou usar só as R1 como base
zcat3 = " ".join([f"{catcmd} {_bn(R1)} | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' > average_length"])
p2 = subprocess.Popen(zcat3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
res_avg = p2.communicate()
if p2.returncode:
    raise Exception(f"system {zcat3} failes: {res_avg[1]}")

# guardar o tamanho médio das reads

average_length2 = 0
# guardar o numero de reads
with open('average_length', "r") as IN:
    for row in IN:
        row = row.rstrip("\n")
        average_length2 = row

gal_file.close()

coverage = (float(average_length2) * soma_reads) / float(genome_size)

tipo_de_resultado = 'coverage'
imprimir = coverage
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

for chau in ('out_R1', 'out_R2', 'average_length'):
    shutil.rmtree(chau)

sys.exit(0)
