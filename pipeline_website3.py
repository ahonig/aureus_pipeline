#!/usr/bin/env python3
# -*- coding: utf8 -*-

import builtins
import functools
import perllib
import re
import sys

_bn = lambda s: '' if s is None else s
_str = lambda s: '' if s is None else str(s)
sys.path[0:0] = ['/opt/pipeline/lib/']
import libs.alinhamento_poli as alinhamento_poli
import libs.alinhamento_poli_enterobacter as alinhamento_poli_enterobacter
import libs.alinhamento_poli_acineto as alinhamento_poli_acineto
import libs.alinhamento_poli_pseudo as alinhamento_poli_pseudo
import libs.alinhamento_outros_pseudo as alinhamento_outros_pseudo
import libs.alinhamento_outros_kleb as alinhamento_outros_kleb
import libs.save_result_mongo as save_result_mongo

perllib.init_package('main')

MLST_result = ''
_d = None
a = ''
b = ''
count2 = 0
docker = perllib.Array()
especieAb = perllib.Array()
especieEc = perllib.Array()
especieKp = perllib.Array()
fasta_outros = ''
fasta_polimixina = ''
identificacao = perllib.Array()
lines = perllib.Array()
lista_acineto = ''
lista_enterobacter = ''
lista_kleb = ''
nearest_sts = ''
preidentificacao = perllib.Array()
resultadoANI = ''
THREADS = "16"  # new one will have 32

sys.argv = dict(perllib.Array(sys.argv)[1:])
builtins.__PACKAGE__ = 'main'
perllib.WARNING = 1
# SKIPPED: use strict;
# SKIPPED: use Scalar::Util qw(looks_like_number);
# use alinhamento_poli_truncation;

# Linha de comando do script perl pipeline_melise_output_gal.pl  <caminho do diretorio Output, com resultados de spades/unicycler e prokka ex:/home/melise/Output_run18-06.02.21> <nome da amostra na pasta Output ex:27563_S12> <caminho onde esta instalado o abricate  ex:./abricate/bin/abricate>
# <arquivo da tabela excel formato .xls onde vão ser impressos os resultados> <diretorio onde esta instalado o kmer-db  ex:/home/melise/kmer-db>  <caminho do diretorio onde esta instalado o mlst ex:/home/melise/identificar_clones> <caminho do DB com as seqs de ptn de resistencia mutacional a
# polimixina ex:/home/melise/resis_poli> <caminho do DB com as seqs de ptn de resistencia mutacional a outros antibioticos ex:/home/melise/outrasMut> <caminho kraken ex:kraken> <caminho unicycle ex:unicycler> <caminho R1> <caminho R2>

# Guardar o caminho do diretorio "Output", com o resultado da montagem, anotacao e blast para todas as amostra. ex /home/melise/Output
caminho1 = sys.argv.get(0)

# Guardar o nome da amostra com o _S*
sample = sys.argv.get(1)
# para pegar o numero da amostra antes do _S1
sample1 = perllib.Array(_str(sample).split('_'))
sample2 = sample1.get(0)
perllib.perl_print(f"Sample: {_bn(sample2)}")

# criar um diretório para a amostra
amostra_dir = perllib.Array(['mkdir', '-p', f"{_bn(caminho1)}/{_bn(sample)}"])
if not (perllib.system(" ".join(amostra_dir)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, amostra_dir))} failes: {perllib.CHILD_ERROR}")

# caminho para onde esta instalado o abricate ex: ./abricate/bin/abricate

caminho_abricate = sys.argv.get(2)

# Caminho para o output onde esta a tabela que sera colocado os resultados
caminho_output = sys.argv.get(3)

# entrar com o caminho da pastar onde esta instalado o kmer-db  /home/melise/kmer-db
kmerdb_install = sys.argv.get(4)

# entrar com o caminho da pastar onde esta instalado o mlst. ex: /home/melise/identificar_clones
mlst_install = sys.argv.get(5)

# Caminho para banco de dados com sequências de genes de R a polimixina
db_polimixina = sys.argv.get(6)

# Caminho para banco de dados com sequências de genes de R mutacionais a outros antibioticos
db_outrosMut = sys.argv.get(7)

# entrar com o caminho da pastar onde esta instalado o kraken2 ex: kraken2
kraken2_install = sys.argv.get(8)

# entrar com o caminho do unicycler
unicycler = sys.argv.get(9)

# criar um diretorio para unicycler
# my @unicycler_dir = ("mkdir","$caminho1/$sample/unicycler");
#		system(@unicycler_dir) == 0
#			or die "system @unicycler_dir failes: $?";

perllib.perl_print('Parametros: ')
perllib.perl_print(f"caminho: {_bn(caminho1)} ")
perllib.perl_print(f"Sample: {_bn(sample)} ")
perllib.perl_print(f"SAmple2: {_bn(sample2)} ")
perllib.perl_print(f"caminho1: {_bn(caminho1)} ")
perllib.perl_print(f"camino abricate: {_bn(caminho_abricate)} ")
perllib.perl_print(f"camino abricate caminho_output: {_bn(caminho_output)} ")
perllib.perl_print(f"camino abricate kmerdb_install: {_bn(kmerdb_install)} ")
perllib.perl_print(f"mlst install: {_bn(mlst_install)}  ")
perllib.perl_print(f"db polimixina: {_bn(db_polimixina)}  ")
perllib.perl_print(f"db outros mut: {_bn(db_outrosMut)}  ")
perllib.perl_print(f"kraken2_install: {_bn(kraken2_install)}  ")
perllib.perl_print(f"unicycler: {_bn(unicycler)} ")

R1 = sys.argv.get(10)
R2 = sys.argv.get(11)

##################################################
# rodar unicycler
if False:
    unicycler_exe = perllib.Array(['unicycler', '-1', f"{_bn(R1)}", '-2', f"{_bn(R2)}", '-o', f"{_bn(caminho1)}/{_bn(sample)}/unicycler", '--min_fasta_length', '500', '--mode', 'conservative', '-t', THREADS, '--spades_path', '/opt/SPAdes-3.13.0-Linux/bin/spades.py'])
    if not (perllib.system(" ".join(unicycler_exe)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, unicycler_exe))} failes: {perllib.CHILD_ERROR}")

# arquivo assembly.fasta da amostra

montagem = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta"

###############################################################################
# rodar prokka
if False:
    prokka_exe = perllib.Array(['prokka', '--outdir', f"{_bn(caminho1)}/{_bn(sample)}/prokka", '--prefix', 'genome', f"{montagem}", '--force', "--cpus", "0"])  # --cpus 0 is ALL
    if not (perllib.system(" ".join(prokka_exe)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, prokka_exe))} failes: {perllib.CHILD_ERROR}")

# EXCLUSIVO DO PIPELINE OUTPUT

# variavel para guardar os tipos de resultados que devem ser indicados
# ex: checkm; especie; contigs; resfinder; VFDB; plasmid; mlst; mutacoes_poli; mutacoes_outra

tipo_de_resultado = None
# o que imprimir desse resultado
imprimir = None

#####################################################################################
# PARA IMPRIMIR O RESULTADO EM UM ARQUIVO TXT PARA GAL

output = 'resultado_gal.txt'
if not ((OUT2 := perllib.open_(_str(output), 'a'))):
    perllib.die('Nao foi possivel abrir output\n')

##printar no arquivo final o nome da amostra na primeira linha

perllib.perl_print(f"\nAmostra {_bn(sample)}\nResultados relevantes do sequenciamento do genoma total (WGS):", file=OUT2)

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
if False:
    copy_contigs = perllib.Array(['mkdir', '-p', 'checkM_bins', '&&', 'cp', f"{montagem}", 'checkM_bins/'])  # creating dir first!
    if not (perllib.system(" ".join(copy_contigs)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, copy_contigs))} failes: {perllib.CHILD_ERROR}")

    # rodar o CheckM

    checkM = perllib.Array(['checkm', 'lineage_wf', '-x', 'fasta', 'checkM_bins', 'checkM_bins', "--threads", THREADS, "--pplacer_threads", THREADS])
    if not (perllib.system(" ".join(checkM)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, checkM))} failes: {perllib.CHILD_ERROR}")

    checkM_qa = perllib.Array(['checkm', 'qa', '-o', '2', '-f', f"checkM_bins/{_bn(sample)}_resultados", '--tab_table', 'checkM_bins/lineage.ms', 'checkM_bins', "--threads", THREADS])
    if not (perllib.system(" ".join(checkM_qa)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, checkM_qa))} failes: {perllib.CHILD_ERROR}")

    # apagar arquivos gerados, deixando apenas resultados

    rm10 = perllib.Array(['rm', '-r', 'checkM_bins/bins'])

    if not (perllib.system(" ".join(rm10)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm10))} failes: {perllib.CHILD_ERROR}")

    rm11 = perllib.Array(['rm', '-r', 'checkM_bins/storage'])

    if not (perllib.system(" ".join(rm11)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm11))} failes: {perllib.CHILD_ERROR}")

    rm12 = perllib.Array(['rm', 'checkM_bins/assembly.fasta'])

    if not (perllib.system(" ".join(rm12)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm12))} failes: {perllib.CHILD_ERROR}")

    rm13 = perllib.Array(['rm', 'checkM_bins/lineage.ms'])

    if not (perllib.system(" ".join(rm13)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm13))} failes: {perllib.CHILD_ERROR}")

    rm14 = perllib.Array(['rm', '-r', 'checkM_bins/checkm.log'])

    if not (perllib.system(" ".join(rm14)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm14))} failes: {perllib.CHILD_ERROR}")

# Ler o arquivo do resultado e imprimir o que interessa na tabela

resultado_checkM = f"checkM_bins/{_bn(sample)}_resultados"
# print "$file\n";
if not ((IN_check := perllib.open_(resultado_checkM, 'r'))):
    perllib.die(f"File {resultado_checkM} not open\n")

# pegar o resultado da contaminacao

contaminacao = 0

perllib.perl_print('Salvando resultado no mongo relatorios')

genome_size = None

cabecera = 'S'
#while (row := perllib.readline_full(IN_check, 'IN_check')):
for row in IN_check:
    # remove \n of the line end
    row = row.rstrip("\n")
    # separar as colunas do arquivo em elementos em um array
    if cabecera == 'N':
        lines = row.split("\t")  # perllib.Array(_str(row).split('\\t'))
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

perllib.perl_print('rodar o kraken', end="")
if False:
    kraken = perllib.Array([f"{_bn(kraken2_install)}/kraken2", '--db', f"{_bn(kraken2_install)}/minikraken2_v2_8GB_201904_UPDATE", '--use-names', '--paired', f"{_bn(R1)}", f"{_bn(R2)}", '--output', 'out_kraken', "--threads", THREADS])
    if not (perllib.system(" ".join(kraken)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, kraken))} failes: {perllib.CHILD_ERROR}")

    out_kraken = 'out_kraken'
    if not ((IN6 := perllib.open_(out_kraken, 'r'))):
        perllib.die(f"File {out_kraken} not open\n")

# I HAVE DOUBTS ####
# ordenado = perllib.Array()
#
# while (row := perllib.readline_full(IN6, 'IN6')):
#     # remove \n of the line end
#     row = row.rstrip("\n")
#     # separar e pegar espécie
#     out_kraken2 = perllib.Array(perllib.split(r'	', _str(row)))
#     # print "$out_kraken2[2]\n";
#     ordenado.extend(perllib.make_list(out_kraken2.get(2)))
#
# ordenado2 = perllib.Array()
#
# for n in ordenado:
#     if ([_m[0] for _m in (re.finditer(re.compile(r'(\w*\s\w*\s).*\(taxid.*\)', re.I), _str(n)))]) or ([_m[0] for _m in (re.finditer(re.compile(r'(\w+\ssp\.).*\(taxid.*\)', re.I), _str(n)))]) or ([_m[0] for _m in (re.finditer(re.compile(r'^(\w+)\s\(taxid.*\)', re.I), _str(n)))]):  # mod 20.05022
#         ordenado2.append(_m.group(1))
# DOUBT ENDS HERE ####

if False:
    ordenado = []
    for row in IN6:
        row = row.rstrip('\n')
        out_kraken2 = row.split(r'	')
        ordenado.append(out_kraken2[2])

    ordenado2 = []
    for n in ordenado:
        a = re.search(r'(\w*\s\w*\s).*\(taxid.*\)', n, re.IGNORECASE)
        b = re.search(r'(\w+\ssp\.).*\(taxid.*\)', n, re.IGNORECASE)
        c = re.search(r'^(\w+)\s\(taxid.*\)', n, re.IGNORECASE)
        for x in (a, b, c):
            if not x:
                continue
            ordenado2.append(x.group(1))

# contar qts vezes aparece cada especie

# I HAVE DOUBTS ####
# count_ordenado2 = perllib.Hash()
# for _d in ordenado2:
#     count_ordenado2[_s0] = perllib.num(count_ordenado2.get(_s0 := _str(_d))) + 1
# # ordenar as repeticoes
#
# repeticoes = perllib.Array(sorted(perllib.Array(count_ordenado2.keys()), key=functools.cmp_to_key(lambda a, b: perllib.spaceship(perllib.num(count_ordenado2.get(_str(a))), perllib.num(count_ordenado2.get(_str(b)))))))
# # repeticoes = sorted(count_ordenado2.keys(), key=lambda x: count_ordenado2[x])
# # pegar as maiores repeticoes
# maior_repeticao = repeticoes.get(-1)
# segunda_repeticao = repeticoes.get(-2)
# DOUBT ENDS HERE ####

if False:
    count_ordenado2 = {}
    for item in ordenado2:
        if item not in count_ordenado2:
            count_ordenado2[item] = 0
        count_ordenado2[item] += 1

    repeticoes = sorted(count_ordenado2.keys(), key=lambda x: count_ordenado2[x])
maior_repeticao = "Pseudomonas aeruginosa"  # repeticoes[-1]
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
resultado_final_especie = ''  # mod 11.05.22
# resultado que sera impresso
printar_especies = ''  # mod 11.05.22
# o que usar para mlst
especie_mlst = ''  # mod 11.05.22

if (re.findall(re.compile(r'\w+\s\w+', re.I), check_especies)):  # mod 11.05.22
    identificar_especie = perllib.Array(check_especies.split())
    genero = f"{_bn(identificar_especie.get(0, ''))}"
    especie = f"{_bn(identificar_especie.get(1, ''))}"

    # print "$genero\n$especie\n";

    # Associar o nome da especie ao banco do mlst e gravar o nome que sera dado como resultado final
    resultado_final_especie = f"{genero}{especie}"
    # $printar_especies = $resultado_final_especie;
    # $especie_mlst = "";

else:
    printar_especies = check_especies
# mod ate aqui 20.05.22

#######################################################################################################################
# Sequencia para verificar mutacoes pontuais usando subrotinas proprias

# guardar o resultado das mutacoes para polimixina

result2 = perllib.Array()
# guardar o resultado das mutacoes para outros antibioticos
result3 = perllib.Array()
# guardar resultados dos fatores de virulencia
vfdb = perllib.Array()

while (resultado_final_especie := f"{genero}{especie}"):
    perllib.perl_print(f"resultado_final_especie: {resultado_final_especie}", end="")
    if resultado_final_especie == 'Pseudomonasaeruginosa':
        especie_mlst = 'paeruginosa'
        printar_especies = 'Pseudomonas aeruginosa'
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_pseudo_poli.fasta"
        result2 = perllib.Array(alinhamento_poli_pseudo.poli_pseudo(montagem, fasta_polimixina, sample, THREADS))
        fasta_outros = f"{_bn(db_outrosMut)}/proteins_outrasMut_pseudo.fasta"
        result3 = perllib.Array(alinhamento_outros_pseudo.outros_pseudo(montagem, fasta_outros, sample, THREADS))
        break
    elif resultado_final_especie == 'Klebsiellapneumoniae':
        especie_mlst = 'kpneumoniae'
        # $printar_especies = "Klebsiella pneumoniae";
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_kleb_poli.fasta"
        result2 = perllib.Array(alinhamento_poli.poli(montagem, fasta_polimixina, sample, THREADS))
        # @result2 = &alinhamento_poli_truncation::poli($montagem,$fasta_polimixina,$sample);
        fasta_outros = f"{_bn(db_outrosMut)}/proteins_outrasMut_kleb.fasta"
        result3 = perllib.Array(alinhamento_outros_kleb.outros_kleb(montagem, fasta_outros, sample, THREADS))

        perllib.perl_print('Para FastANI', end="")
        # Abrir o arquivo lista
        lista = '/opt/genomas_enterobacter/kleb_database/lista-kleb'  # CAMBIAR

        fastani = perllib.Array(['/opt/FastANI/fastANI', '-q', f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", '--rl', f"{lista}", '-o', f"{_bn(sample)}_out-fastANI", "--threads", THREADS])
        if not (perllib.system(" ".join(fastani)) == 0):
            perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, fastani))} failes: {perllib.CHILD_ERROR}")

        # abrir output
        # Abrir o arquivo do output de distancia
        resultadoANI = f"{_bn(sample)}_out-fastANI"
        if not ((IN7 := perllib.open_(resultadoANI, 'r'))):
            perllib.die(f"File {resultadoANI} not open\n")

        # array para guardar especies

        especieKp = perllib.Array()

        count2 = 0

        perllib.perl_print('resultado do fastANI', end="")
        # Ler o arquivo distancia
        # while (_d := perllib.readline_full(IN7, 'IN7')):
        for _d in IN7:
            perllib.perl_print(f"{_bn(_d)}")
            # remove \n of the line end
            count2 += 1
            _d = _d.rstrip("\n")
            if (perllib.num((count2 := count2 + 1) - 1) == 1):
                especieKp = perllib.Array(_str(_d).split('\\t'))

        preidentificacao = perllib.Array(perllib.split(r'/', _str(especieKp.get(1))))
        identificacao = perllib.Array(perllib.split(r'\.', _str(preidentificacao.get(-1))))
        printar_especies = identificacao.get(0)
        break
    elif resultado_final_especie == 'Escherichiacoli':
        especie_mlst = 'ecoli'
        printar_especies = 'Escherichia coli'
        break
    elif resultado_final_especie == 'Staphylococcusaureus':
        especie_mlst = 'saureus'
        printar_especies = 'Staphylococcus aureus'
        break
    elif resultado_final_especie == 'Pseudomonasputida':
        especie_mlst = 'pputida'
        printar_especies = 'Pseudomonas putida'
        break
    elif resultado_final_especie == 'Listeriamonocytogenes':  # modificado 19.11.21
        especie_mlst = 'lmonocytogenes'
        printar_especies = 'Listeria monocytogenes'
        break
    elif resultado_final_especie == 'Acinetobacterbaumannii':
        especie_mlst = 'abaumannii_2'
        # $printar_especies = "Acinetobacter baumannii";
        fasta_polimixina = f"{_bn(db_polimixina)}/proteins_acineto_poli.fasta"
        result2 = perllib.Array(alinhamento_poli_acineto.poli_acineto(montagem, fasta_polimixina, sample, THREADS))

        perllib.perl_print('Para FastANI', end="")
        # Abrir o arquivo lista
        lista = '/opt/genomas_enterobacter/fastANI_acineto/list-acineto'  # CAMBIAR

        fastani = perllib.Array(['/opt/FastANI/fastANI', '-q', f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", '--rl', f"{lista}", '-o', f"{_bn(sample)}_out-fastANI", "--threads", THREADS])
        if not (perllib.system(" ".join(fastani)) == 0):
            perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, fastani))} failes: {perllib.CHILD_ERROR}")

        # abrir output
        # Abrir o arquivo do output de distancia
        resultadoANI = f"{_bn(sample)}_out-fastANI"
        if not ((IN7 := perllib.open_(resultadoANI, 'r'))):
            perllib.die(f"File {resultadoANI} not open\n")

        # array para guardar especies

        especieAb = perllib.Array()

        count2 = 0

        # Ler o arquivo distancia
        #while (_d := perllib.readline_full(IN7, 'IN7')):
        for _d in IN7:
            # remove \n of the line end
            count2 += 1
            _d = _d.rstrip("\n")
            if (perllib.num((count2 := count2 + 1) - 1) == 1):
                especieAb = perllib.Array(_str(_d).split('\\t'))

        preidentificacao = perllib.Array(perllib.split(r'/', _str(especieAb.get(1))))
        identificacao = perllib.Array(perllib.split(r'\.', _str(preidentificacao.get(-1))))
        printar_especies = identificacao.get(0)
        break
    elif resultado_final_especie == 'Enterococcusfaecalis':
        especie_mlst = 'efaecalis'
        printar_especies = 'Enterococcus faecalis'
        break
    elif resultado_final_especie == 'Klebsiellaoxytoca':
        especie_mlst = 'koxytoca'
        printar_especies = 'Klebsiella oxytoca'
        break
    elif resultado_final_especie == 'Enterococcusfaecium':
        especie_mlst = 'efaecium'
        printar_especies = 'Enterococcus faecium'
        break
    elif resultado_final_especie == 'Serratiamarcescens':
        especie_mlst = 'Nao disponivel'
        printar_especies = 'Serratia marcescens'
        break
    elif resultado_final_especie == 'Providenciastuartii':
        especie_mlst = 'Nao disponivel'
        printar_especies = 'Providencia stuartii'
        break
    elif resultado_final_especie.lower() in ("enterobactercloacae", "enterobacterhormaechei", "enterobacterasburiae", "enterobacterkobei", "enterobacterroggenkampii", "enterobacterludwigii"):

        perllib.perl_print('Rodar fastANI para subespecie', end="")
        # Abrir o arquivo lista
        lista = '/opt/genomas_enterobacter/fastANI/list_entero'  # CAMBIAR

        fastani = perllib.Array(['/opt/FastANI/fastANI', '-q', f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", '--rl', f"{lista}", '-o', f"{_bn(sample)}_out-fastANI", "--threads", THREADS])
        if not (perllib.system(" ".join(fastani)) == 0):
            perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, fastani))} failes: {perllib.CHILD_ERROR}")

        # abrir output
        # Abrir o arquivo do output de distancia do kmer-db
        resultadoANI = f"{_bn(sample)}_out-fastANI"
        if not ((IN7 := perllib.open_(resultadoANI, 'r'))):
            perllib.die(f"File {resultadoANI} not open\n")

        # array para guardar especies

        especieEc = perllib.Array()
        count2 = 0

        perllib.perl_print('EspecieEC', end="")

        # Ler o arquivo distancia
        #while (_d := perllib.readline_full(IN7, 'IN7')):
        for _d in IN7:
            # remove \n of the line end
            count2 += 1
            _d = _d.rstrip("\n")
            if (perllib.num((count2 := count2 + 1) - 1) == 1):
                especieEc = perllib.Array(_str(_d).split('\\t'))

        especie_mlst = 'ecloacae'
        preidentificacao = perllib.Array(perllib.split(r'/', _str(especieEc.get(1))))
        identificacao = perllib.Array(perllib.split(r'\.', _str(preidentificacao.get(-1))))
        printar_especies = identificacao.get(0)

        if (re.findall(re.compile(r'Enterobacter_cloacae_subsp_cloacae', re.I), _str(printar_especies))):
            fasta_polimixina = f"{_bn(db_polimixina)}/proteins_Ecloacae_poli.fasta"
            result2 = perllib.Array(alinhamento_poli_enterobacter.poli_enterobacter(montagem, fasta_polimixina, sample, THREADS))

        break
        # imprimir a espécie
        # print "$especieEc[$index]\n";
    else:
        # mod 20.05.22
        printar_especies = f"{genero} {especie}"  # mod 10.05.22
        especie_mlst = 'Nao disponivel'  # mod 26.08.22
        break  # mod 20.05.22
    # mod 20.05.22

# print "$especie_mlst ";

perllib.perl_print('contaminacao...', end="")
# printar no arquivo final o nome da especie
if perllib.num(contaminacao) <= 10:
    # print OUT2 "$printar_especies\t";
    tipo_de_resultado = 'especie'
    imprimir = printar_especies
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    perllib.perl_print(f"Espécie identificada: {_bn(imprimir)}", file=OUT2)

if perllib.num(contaminacao) > 10:
    # print OUT2 "$printar_especies\t";
    tipo_de_resultado = 'especie'
    imprimir = f"{maior_repeticao} {_bn(count_ordenado2.get(maior_repeticao, ''))} {segunda_repeticao} {_bn(count_ordenado2.get(segunda_repeticao, ''))}"
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    perllib.perl_print(f"Espécie: CONTAMINAÇÃO {_bn(imprimir)}", file=OUT2)

# else {
#        	print OUT2 "$maior_repeticao $count_ordenado2{$maior_repeticao} $segunda_repeticao $count_ordenado2{$segunda_repeticao}\t";
# }

###############################################################################################
# Rodar ABRICATE
# Para resistencia usando o ResFinder (porque so tem resistencia adquirida)
abricante_exe = perllib.Array([f"{_bn(caminho_abricate)}", "--db", "resfinder", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", f"{_bn(sample)}_outAbricateRes", '--threads', THREADS])
if not (perllib.system(" ".join(abricante_exe)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, abricante_exe))} failes: {perllib.CHILD_ERROR}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# abrir o resultado do abricate_resfinder Contigs
result_abricate_resfinder = f"{_bn(sample)}_outAbricateRes"
if not ((IN8 := perllib.open_(result_abricate_resfinder, 'r'))):
    perllib.die(f"File {result_abricate_resfinder} not open\n")

# criar um array para guardar so os que tiverem identidade e cobertura alta

selected = perllib.Array()
rowQty = 0
#while (row := perllib.readline_full(IN8, 'IN8')):
for row in IN8:
    rowQty += 1
    # remove \n of the line end
    row = row.rstrip("\n")
    # print "$row\n";
    # separar as colunas do arquivo em elementos em um array
    lines = row.split("\t")  # perllib.Array(_str(row).split('\\t'))
    # print "$lines[10]\n";
    # printar no array selected so quem tem identidade maior que 90
    if rowQty > 4:
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected.append(f"{_bn(row)}")
            perllib.perl_print(f"{_bn(row)}")
        # if para o operon de Van no Resfinder

        if re.match(r'Van.*', lines[5], re.I):
            selected.append(f"{_bn(row)}")
            perllib.perl_print(f"{_bn(row)}")

# Analisar o array dos contigs
# ler o array selected


select_imprimir = perllib.Array()
out_blast = ''

# criar um @ para cada classe de antibioticos
genes = perllib.Array()

# print gal
perllib.perl_print('Genes de resistência encontrados no WGS: ', file=OUT2)

for n_l in selected:
    # print "$n\n";
    # separar as colunas do arquivo em elementos de um array
    lines_blast = n_l.split("\t")  # perllib.Array(_str(n_l).split('\\t'))
    # concatenar os resultado
    out_blast = lines_blast[5] + ' ' + '(' + 'ID' + ':' + lines_blast[10] + ' ' + 'COV_Q:' + lines_blast[9] + ' ' + 'COV_DB:' + lines_blast[6] + ')'
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
abricante_exe = perllib.Array([f"{_bn(caminho_abricate)}", "--db", "vfdb", f"{_bn(caminho1)}/{_bn(sample)}/prokka/genome.ffn", ">", f"{_bn(sample)}_outAbricateVFDB", '--threads', THREADS])
if not (perllib.system(" ".join(abricante_exe)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, abricante_exe))} failes: {perllib.CHILD_ERROR}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# abrir o resultado do abricate_VFDB
result_abricate_VFDB = f"{_bn(sample)}_outAbricateVFDB"
if not ((IN10 := perllib.open_(result_abricate_VFDB, 'r'))):
    perllib.die(f"File {result_abricate_VFDB} not open\n")

# criar um array para guardar so os que tiverem identidade e cobertura alta

selected3 = perllib.Array()

select_imprimir3 = perllib.Array()
out_blast3 = ''
qty = 0
#while (row := perllib.readline_full(IN10, 'IN10')):
for row in IN10:
    qty += 1
    # remove \n of the line end
    row = row.rstrip("\n")
    # print "$row\n";
    # separar as colunas do arquivo em elementos em um array
    lines = row.split("\t")  # perllib.Array(_str(row).split('\\t'))
    # print "$lines[10]\n";
    # printar no array selected so quem tem coverage maior que 90
    if qty > 4:
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected3.append(f"{_bn(row)}")

# ler o array selected

for n_l in selected3:
    # print "$n\n";
    # separar as colunas do arquivo em elementos de um array
    lines_blast = n_l.split("\t")  # perllib.Array(_str(n_l).split('\\t'))
    # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
    out_blast3 = lines_blast[1] + ':' + ' ' + lines_blast[5] + ' ' + lines_blast[13] + ' ' + 'ID' + ':' + lines_blast[10] + ' ' + 'COV_Q:' + lines_blast[9] + ' ' + 'COV_DB:' + lines_blast[6] + '|' + ' '
    select_imprimir3.append(out_blast3)

# zerar a variarvel para concatenar
imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in select_imprimir3:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'VFDB'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

rm17 = perllib.Array(['rm', f"{_bn(sample)}_outAbricateVFDB"])

if not (perllib.system(" ".join(rm17)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm17))} failes: {perllib.CHILD_ERROR}")

#########################################################################################################
# Rodar abricate para PlasmidFinder
abricante_exe = perllib.Array([f"{_bn(caminho_abricate)}", "--db", "plasmidfinder", f"{_bn(caminho1)}/{_bn(sample)}/unicycler/assembly.fasta", ">", f"{_bn(sample)}_outAbricatePlasmid", '--threads', THREADS])
if not (perllib.system(" ".join(abricante_exe)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, abricante_exe))} failes: {perllib.CHILD_ERROR}")
# system("$caminho_abricate --db resfinder $caminho1/$sample/ContigsMenor/prokka/genome.ffn > $sample\_outAbricateResMenor");

# abrir o resultado do abricate_Plasmid
result_abricate_Plasmid = f"{_bn(sample)}_outAbricatePlasmid"
if not ((IN11 := perllib.open_(result_abricate_Plasmid, 'r'))):
    perllib.die(f"File {result_abricate_Plasmid} not open\n")

# criar um array para guardar so os que tiverem identidade e cobertura alta

selected4 = perllib.Array()

select_imprimir4 = perllib.Array()
out_blast4 = None
qxty = 0
#while (row := perllib.readline_full(IN11, 'IN11')):
for row in IN11:
    qxty += 1
    # remove \n of the line end
    row = row.rstrip("\n")
    # print "$row\n";
    # separar as colunas do arquivo em elementos em um array
    lines = row.split("\t")  # perllib.Array(_str(row).split('\\t'))
    # print "$lines[10]\n";
    # printar no array selected so quem tem coverage maior que 90
    if qxty > 4:
        if (float(lines[9]) > 90.0) and (float(lines[10]) > 90.0):
            selected4.append(f"{_bn(row)}")
# ler o array selected

if len(selected4) == 0:
    # print OUT2 "Nao encontrado\t";
    tipo_de_resultado = 'plasmid'
    imprimir = 'Not found'
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    perllib.perl_print(f"Plasmídeos encontrados: {_bn(imprimir)}", file=OUT2)
else:
    perllib.perl_print('Plasmídeos encontrados:', file=OUT2)
    for n_l in selected4:
        # print "$n\n";
        # separar as colunas do arquivo em elementos de um array
        lines_blast = n_l.split("\t")  # perllib.Array(_str(n_l).split('\\t'))
        # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
        out_blast4 = lines_blast[5] + 'ID' + ':' + lines_blast[10] + ' ' + 'COV_Q:' + lines_blast[9] + ' ' + 'COV_DB:' + lines_blast[6] + '|' + ' '
        select_imprimir4.append(out_blast4)
        perllib.perl_print(f"{_bn(lines_blast[5])}", file=OUT2)
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

rm18 = perllib.Array(['rm', f"{_bn(sample)}_outAbricatePlasmid"])

if not (perllib.system(" ".join(rm18)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm18))} failes: {perllib.CHILD_ERROR}")

################################################################################################

perllib.perl_print(f"Rodar o MLST {especie_mlst}", end="")

# se nao tem mlst disponivel, ai tem que avisar
if (especie_mlst == 'Nao disponivel') or (especie_mlst == ''):  # mod 26-08-22
    # print OUT2 "Nao disponivel\t";
    tipo_de_resultado = 'mlst'
    imprimir = 'Not available for this species'  # mod 26.08.22
    save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
    # para o gal
    perllib.perl_print(f"Clone ST {_bn(imprimir)} (determinado por MLST)", file=OUT2)
else:
    # mod 26-08-22
    docker = perllib.Array(['docker', 'run', '--rm', '-i', '-v', f"{_bn(mlst_install)}/mlst_db:/database", '-v', f"{_bn(caminho1)}/{_bn(sample)}/unicycler:/workdir", 'mlst', '-i', 'assembly.fasta', '-o', '.', '-s', f"{especie_mlst}"])

    # rodar o mlst
    if not (perllib.system(" ".join(docker)) == 0):
        perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, docker))} failes: {perllib.CHILD_ERROR}")

    # abrir o resultado do mlst

    MLST_result = f"{_bn(caminho1)}/{_bn(sample)}/unicycler/data.json"
    if not ((IN3 := perllib.open_(MLST_result, 'r'))):
        perllib.die(f"File {MLST_result} not open\n")

# mod

ST = None
perllib.perl_print('ler o resultado do mlst', end="")
#while (row3 := perllib.readline_full(IN3, 'IN3')):
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
        perllib.perl_print(f"Clone ST {_bn(imprimir)} (determinado por MLST)", file=OUT2)
    m = re.match(r'nearest_sts":\s"((\d*,)*\d*)".*', row3, re.IGNORECASE)
    if m:
        nearest_sts = m.group(1)
        # print OUT2 "Nearest $nearest_sts\t";
        tipo_de_resultado = 'mlst'
        imprimir = f"Nearest {nearest_sts}"
        save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
        # para o gal
        perllib.perl_print(f"Clone ST {_bn(imprimir)} (determinado por MLST)", file=OUT2)
    m = re.match(r'.*sequence_type":\s"(Unknown)".*', row3, re.IGNORECASE)
    if m:
        ST = m.group(1)
        # print OUT2 "Unknown\t";
        tipo_de_resultado = 'mlst'
        imprimir = 'Unknown'
        save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
        # para o gal
        perllib.perl_print(f"Clone ST {_bn(imprimir)} (determinado por MLST)", file=OUT2)

# zerar a variarvel para concatenar

imprimir = ''

# colocar o resultado que estava salvo no @select_imprimir  na variável $imprimir
for _d in result2:
    imprimir = _str(imprimir) + _str(_d)
    imprimir = _str(imprimir) + '<br>'

tipo_de_resultado = 'mutacoes_poli'
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)
# para o gal
perllib.perl_print(f"Mutações polimixina: {_bn(imprimir)}\n", file=OUT2)

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
perllib.perl_print('rodar coverage', end="")
zcat = perllib.Array([f"echo $(zcat {_bn(R1)} | wc -l)/4 | bc > out_R1 "])
# print "@zcat\n";
if not (perllib.system(" ".join(zcat)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, zcat))} failes: {perllib.CHILD_ERROR}")

out_R1 = 'out_R1'
if not ((IN := perllib.open_(out_R1, 'r'))):
    perllib.die(f"File {out_R1} not open\n")

n_reads1 = None
# guardar o numero de reads
#while (row := perllib.readline_full(IN, 'IN')):
for row in IN:
    row = row.rstrip("\n")
    n_reads1 = row

# o mesmo para o arquivo R2

zcat2 = perllib.Array([f"echo $(zcat {_bn(R2)} | wc -l)/4 | bc > out_R2"])
if not (perllib.system(" ".join(zcat2)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, zcat2))} failes: {perllib.CHILD_ERROR}")

out_R2 = 'out_R2'
if not ((IN2 := perllib.open_(out_R2, 'r'))):
    perllib.die(f"File {out_R2} not open\n")

n_reads2 = None
# guardar o numero de reads
#while (row2 := perllib.readline_full(IN2, 'IN2')):
for row2 in IN2:
    row2 = row2.rstrip("\n")
    n_reads2 = row2

soma_reads = (perllib.num(n_reads1) + perllib.num(n_reads2))

###calcular tamanho medio das reads, vou usar só as R1 como base
zcat3 = perllib.Array([f"zcat {_bn(R1)} | awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' > average_lenght"])
if not (perllib.system(" ".join(zcat3)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, zcat2))} failes: {perllib.CHILD_ERROR}")

average_length = 'average_lenght'
if not ((IN3 := perllib.open_(average_length, 'r'))):
    perllib.die(f"File {average_length} not open\n")

# guardar o tamanho médio das reads

average_length2 = 0
#while (row3 := perllib.readline_full(IN3, 'IN3')):
for row3 in IN3:
    row3 = row3.rstrip("\n")
    average_length2 = row3

perllib.close_(IN)
perllib.close_(IN2)
perllib.close_(IN3)

coverage = (perllib.num(average_length2) * soma_reads) / perllib.num(genome_size)

tipo_de_resultado = 'coverage'
imprimir = coverage
save_result_mongo.save_result(caminho_output, sample2, tipo_de_resultado, imprimir)

rm19 = perllib.Array(['rm', 'out_R1'])

if not (perllib.system(" ".join(rm19)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm19))} failes: {perllib.CHILD_ERROR}")

rm20 = perllib.Array(['rm', 'out_R2'])

if not (perllib.system(" ".join(rm20)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm20))} failes: {perllib.CHILD_ERROR}")

rm21 = perllib.Array(['rm', 'average_lenght'])

if not (perllib.system(" ".join(rm21)) == 0):
    perllib.die(f"system {perllib.LIST_SEPARATOR.join(map(_str, rm21))} failes: {perllib.CHILD_ERROR}")

sys.exit()
