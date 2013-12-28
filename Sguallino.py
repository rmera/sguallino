#!/usr/bin/python
#############################SGUALLINO 2.0-b2, o SuperGuallino#####################################
#To the long life of the Ven. Khenpo Puntzok Tenzin Rinpoche

 # This software is Copyright (C) 2005 by Raul Mera-Adasme, rmeraa+soft@gmail.com
 # 
 #This program is free software; you can redistribute it and/or
 #modify it under the terms of the GNU General Public License
 #as published by the Free Software Foundation; either version 2
 #of the License, or (at your option) any later version. 

 #This program is distributed in the hope that it will be useful,
 #but WITHOUT ANY WARRANTY; without even the implied warranty of
 #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #GNU General Public License for more details.

 #You should have received a copy of the GNU General Public License
 #along with this program; if not, write to the Free Software
 #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



#Sorry, Spanish.
#I am not very likely to ever translate this.

########################################################################
#esta es una reescritura en python del programa GUALLINO, incorporando otras
#caracteristicas que lo convierten en un completo evaluador de alignments. (estas estan listadas mas arriba)
#actualmetnte:
#pondera el alignment segun la semejanza en polaridad tamano y reactividad de los dos aminoacidos alineados
#superpone al alignment la estructura secundaria experimental para una proteina template y una 
##estructura secundaria teorica para la proteina a modelar para ver en que parte de la estructura quedan los gap
#pondera la ocurrencia en proteinas de cada aminoacido de la proteína a modelar.
###este programa recibe un archivo de alignment en formato fasta, una linea por 
#secuencia, la primera secuencia es la problem y la segunda es la que se ha usado para 
#obtener la estructura secundaria experimental, (esta estructura debe incluirse en un archivo llamado 
#dssp.sec. También se debe incluir una prediccion de estructura secundaria para la primera 
#secuencia en un archivo llamado predic.sec)
#Solo las funciones con una P antes de su nombre, como P_filter son publicas. 
#(hay una clase en la que esto ultimo no se cumple ARREGLAR)
#Produce un PDB que indica la calidad del alignment, ademas de los lugares donde hay insreciones y 
#deleciones (ver la documentación del programa, o los comentarios de abajo para más detalles)
#UPDATE: SE HA CAMBIADO EL ORDEN QUE DEBEN LLEVAR LAS SECUENCIAS EN EL ALINEAMIENTO:
#AHORA LA SECUENCIA TEMPLADO IRÁ EN PRIMER LUGAR Y LA PROBLEMA EN SEGUNDO LUGAR EN EL ARCHIVO
#ALINEAMIENTO.fasta
###########################################


#USAGE: ./Sguallino.py [-p] [-pdb file_in.pdb] [-pedantic]  [-pdbout file_out.pdb] [-sequence] [-fileout file_out]


###Opcion -p dispays the problem sequence numering instead of the template.
###Opcion -fileout followed by a filename gives the chosen name to the output, which by default is called OUT. 
###Option -matrix file_matrix allows to specify a file containing a comparison matrix. 



#en cada linea (el archivo debe tener tres) debe tener una clasificacion de aminoácidos 
###La opción -pdb seguida por un nombre de archivo hace que se genere a partir de ese archivo un 
#pdb con factores de temperatura que indican el ajuste de los aa de la secuencia problem y templado
#se diferencian por defecto: sin correspondencia, correspondencia de tamaño, correspondencia de polaridad
#correspondencia de tamaño y polaridad, en orden decreciente de b_factor estos se pueden observar con PyMOL
#y los peores ajustes se verán en rojo. Además este PDB incluye atomos de K+ que marcan sitios de deleciones
#y atomos de I- que marcan inserciones (desde el punto de vista de la secuencia problem).
#todas las opciones siguientes requieren que esta esté activada.
###La opción pedantic agrega a las distinciones de alignment del PDB producido por la opcion -pdb
#la de ajuste perfecto, que se verá como azúl oscuro en pymol.
###La opcion -pdbout seguida de un nombre de archivo da un nombre al pdb de salida de la opción -pdb.
#sin esta opcion, el pdb de salida se llamara alin_pdb.pdb
###Con la opción -frequences los identificadores de cadena de cada aminoácido se reemplacen por un indicador 
#de su abundancia relativa en proteínas. Esto permite identificar los aminoácidos mas o menos comunes
#pero impide su correcta visualizacion en los modos "cartoons" y "ribbons" de pymol. Una buena idea
#es correr Sguallino 2 veces y producir 2 PDB diferentes con y sin la opcion
#para verlos superpuestos en PyMOL (basta abrir uno despues de otro). Esta opción deshabilita los indicadores
#de inserciones y deleciones.

import sys
from ShowPdb import alin_pdb
from matrix_selector import matrix_selector
from copy import copy
class filehandle:
	def __init__(self):
		self.Falin=open("./ALIN.fasta")
		self.Fsece=open("./dssp.sec")
		self.Fsect=open("./prediction.sec")
		self.Malin=self.Falin.readlines()
		self.Msece=self.Fsece.readline()
		self.Msect=self.Fsect.readline()
		self.Falin.close()
		self.Fsece.close()
		self.Fsect.close()	

class abund_aa:
	def __init__(self,secuencia):
		self.template=[]
		self.template=list(secuencia)
		#My own Q&D abundance classification: Larger than 6.0, larger than 4.0, larger than 2.0 and the rest
		self.pond_goyira=[["G","A","V","L","S","E"],["I","P","T","N","Q","K","R","D"],["M","Y","H","F"],["W","C"]]
		self.pond=[["A","G","V","L","S"],["I","T","K","D","E"],["R","P","F","N","Q"],["Y","H","C","M","W"]]
		self.out=""
		self.cont=0
		for i in self.template:
			if i=="-":
				self.out=self.out+"-"
			elif self.pond[0].count(i):
				self.out=self.out+"F"
			elif self.pond[1].count(i):
				self.out=self.out+"C"
			elif self.pond[2].count(i):
				self.out=self.out+"U"
			elif self.pond[3].count(i):
				self.out=self.out+"R"
			else:
				self.out=self.out+"X"
				print "error in class abund_aa!!!!!!!!!!!! aa numer:", self.cont
			self.cont=self.cont+1
				

class numero_aa:
#Creates a vector with the number of each aa in the template
	def __init__(self, secuencia):
		self.template=[]
		self.template=list(secuencia)
		self.out=""
		self.cont=0
		self.cont_dec=1
		self.signal=0 #How many digits above 1 had the last-written number.
		for i in self.template:
			
			if i=="-":
				if self.signal>0:
					self.signal=self.signal-1
					continue
				else:
					self.out=self.out+" "
			else:
				self.cont=self.cont+1
				if self.signal>=1:
					self.signal=self.signal-1
					continue
				if self.cont-self.cont_dec==10:
					self.cont_dec=self.cont_dec+10
					self.out=self.out+str(self.cont)
					if self.cont>9 and self.cont<=99:
						self.signal=1
					elif self.cont>99 and self.cont<=999:
						self.signal=2
					elif self.cont>999 and self.cont<=9999:
						self.signal=3
					elif self.cont>9999: #This should be unreachable.
						self.signal=4

				else:
					self.out=self.out+" "
					
									
		

#This code is in Spanish. I was not too good at english with i wrote it, sorry.
#I have tried to change some names and translate the comments.
class sorter:
	def __init__(self, alignment, este, estt):
		self.alignment=alignment
		self.este=este
		self.estt=estt
		self.Malin=[]
	def orden(self,template,struct):
		self.template=template
		self.struct_in=list(struct)
		self.struct_out=[]
		self.salto=0
		self.vueltas=0
		for j in self.template:
			if j=="-":  
				self.struct_out.append(j)
				self.salto=self.salto+1
				self.vueltas=self.vueltas+1
			else:
				self.indice=self.vueltas-self.salto
				if len(self.struct_in)<=self.indice:
					self.struct_out.append("X")
					self.vueltas=self.vueltas+1
					continue
				self.struct_out.append(self.struct_in[self.indice])
				self.vueltas=self.vueltas+1
		if self.struct_out.count("\n")>0:
			self.struct_out.remove("\n")
		return self.struct_out
	def P_filter(self):
		self.nombres_sec=[]
		self.status=0
		self.revisadas=0
		for i in self.alignment:
			self.lista=list(i)
			if self.lista[0]==">":
				self.status=self.status+1
				self.nombres_sec.append(i.replace("\n",""))
				continue
			if self.status<=self.revisadas:
				continue
			self.revisadas=self.revisadas+1
			if self.revisadas==1:
				self.este=self.orden(self.lista,self.este  )
				self.Malin.append(i.replace("\n",""))
			elif self.revisadas==2:
				self.estt=self.orden(self.lista,self.estt)
				self.Malin.append(i.replace("\n",""))
			else:
				self.Malin.append(i.replace("\n",""))
		self.lista=[]

#3 criteria are evaluated: polarity likeness, and secondary structure likeness.
#Met has been included in the non-reactive group, about which i am not completely sure.
#the final number for the criteria is the sum of the ones that apply for a given pair:
#polarity equals 4, beta-sheet 2 and alpha-helix 1. If a pair matches in all criteria the final
#value would be 7.
class evaluador:
	def __init__(self,M,file_matrix):
		self.M=M
		self.Mlistas=[]
		self.identity=[]
		if file_matrix:
			self.lista=matrix_selector(file_matrix)
		else:
			#Polarity size and chemical group.
			#self.lista_ant=[[["G","A","V","L","I","M","P","F","W"],["S","T","N","Q","Y","C"],["K","R","H"],["D","E"]],[["G"],["A"],["V","T","S","D","C","P"],["L","I","N","E"],["M","Q"],["H","F","W","Y","K","R"]],[["P"],["G","A","V","L","I","M"],["F","W","Y",],["C","H","D","G","E"],["S","T"],["N","Q"],["K","R"],]]
			#polaridad, alfa_helix, beta sheet
			self.lista=[[["G","A","V","L","I","M","P","F","W"],["S","T","N","Q","Y","C"],["K","R","H"],["D","E"]],[["A","Q","E","I","L","K","M","F","W","V"],["R","D","C","H","S","T"],["N","G","P","V"]],[["C","Q","H","I","L","M","F","T","W","Y","V"],["A","R","N"],["D","E","G","K","M","P","S"]]]
		for i in range(len(self.M)):
			self.Mlistas.append(list(self.M[i]))
	def compara(self):
		self.sum=0
		self.Msum=[]
		self.MsumT=[]
		for i in range(len(self.Mlistas)):
			if i==0:       #####not proper.
				continue
			self.MsumT=[]
			for j in range(len(self.Mlistas[i])):
				for k in range(len(self.lista)):
					for f in range(len(self.lista[k])):
						self.match=0
						if self.Mlistas[0][j]==self.Mlistas[i][j] and self.Mlistas[i][j] != "-": ###potential issue here?
							self.sum=7.2
							continue
						if k==1:
							if self.lista[k][f].count(self.Mlistas[0][j]) and self.lista[k][f].count(self.Mlistas[i][j]):
								self.match=1
								break
							elif f>0:
								if self.lista[k][f].count(self.Mlistas[0][j]) and self.lista[k][f-1].count(self.Mlistas[i][j]):
									self.match=1
									break
							elif f<5:
								if self.lista[k][f].count(self.Mlistas[0][j]) and self.lista[k][f+1].count(self.Mlistas[i][j]):
									self.match=1
									break
						else:
							if self.lista[k][f].count(self.Mlistas[0][j]) and self.lista[k][f].count(self.Mlistas[i][j]):
								self.match=1
								break
					if self.match==0:
						continue
					if k==0:
						self.sum=self.sum+4
					if k==1:
						self.sum=self.sum+1
					if k==2:
						self.sum=self.sum+2
				self.MsumT.append(self.sum)
				self.sum=0
			self.Msum.append(self.MsumT)	
			#Here I calculate the percentage of identity beteen the problem sequence and the template, using the problem as a reference			
			self.identity.append(100*(float(self.MsumT.count(7.2))/(len(self.Mlistas[i])-self.Mlistas[i].count("-"))))
			print "Total", len(self.Mlistas[i])-self.Mlistas[i].count("-") , "Equals:", self.MsumT.count(7.2), "Id", 100*(float(self.MsumT.count(7.2))/(len(self.Mlistas[i])-self.Mlistas[i].count("-")))  #############################

		
class troceador:
	def __init__(self):
		self.lista=[]
		self.largo=0
		self.t=0
		self.matriz=[]
		self.inf=0
		self.sup=0
	def troc(self):
		self.temp=""
		if (self.largo<=self.sup):
			for i in range(self.inf,self.largo):
				self.temp=self.temp+str(self.lista[i])
			self.matriz.append(self.temp)

		else:
			for i in range(self.inf,self.sup):
				self.temp=self.temp+str(self.lista[i])
			self.matriz.append(self.temp)
	def P_go(self,lista,tamano):
		self.matriz=[]
		self.lista=[]
		self.lista=lista
		self.largo=len(self.lista)
		print self.largo
		self.t=tamano
		self.inf=0
		self.sup=self.t
		for i in range(0,self.largo+80,self.t): #The previous version used to have a self.Inf where the 0 is now.
			self.sup=i+self.t
			self.troc()
			self.inf=self.sup


		
		
					

F=filehandle()
aa=len(F.Malin[1])
print F.Malin[0]
print aa
O=sorter(F.Malin,F.Msece,F.Msect)
O.P_filter()
este=""
estt=""
Frec=abund_aa(O.Malin[0])

if(len(sys.argv)>1 and sys.argv[1]=="-p"): #  Here we decide whether to prinit the prolem or the template sequence.
	num_aa=numero_aa(O.Malin[1])
	text_sec="problem"
else:
	num_aa=numero_aa(O.Malin[0])
	text_sec="template"

if "-matrix" in sys.argv:
	index=1+sys.argv.index("-matrix")
	matrix=sys.argv[index]
else:
	matrix=False
E=evaluador(O.Malin,matrix)
E.compara()

for j in range(len(O.Malin)):
	O.Malin[j]=list(O.Malin[j])

sum=copy(E.Msum[0])###CHANGED # This keeps the decimals lost for E.Msum in the next step.
for j in range(len(E.Msum)):
	E.Msum[j]=E.Msum[j] ##CHANFED IT WAS list(E.Msum[j]) same line 333
	for k in range(len(E.Msum[j])):
		E.Msum[j][k]=int(E.Msum[j][k])

# Here we put together the PDB showing the alignment quality, together with insertions and deletions.
if "-pdb" in sys.argv:
	index=1+sys.argv.index("-pdb")
	pdb_in=sys.argv[index]
	pdb=alin_pdb(pdb_in)
	pdb.leer_snap()
	if "-pdbout" in sys.argv:
		index=1+sys.argv.index("-pdbout")
		pdb_out=sys.argv[index]
	else:
		pdb_out="alin_pdb.pdb"
	if not "-pedantic" in sys.argv:
		sum=E.Msum[0]  # If the "-pedantic" option is not used, decimals are lost and a perfect match is not separated from a match of size and polarity   
	if "-frequences" in sys.argv:
		pdb.do_color(sum,O.Malin[0],list(Frec.out)) # The frequency PDB doesnt have insertions or deletions.
	else:
		pdb.do_color(sum,O.Malin[0],False)
		pdb.add_deletions(O.Malin[0]) #Insertions/deletions are seen from the problem sequence's point of view.
		pdb.add_insertions(O.Malin[0],O.Malin[1])
	pdb.write_alin_pdb(pdb_out)
#END PDB part

for j in range(len(E.Msum)):
	for k in range(len(E.Msum[j])):
		E.Msum[j][k]=int(E.Msum[j][k])
		
	


Malinprint=[]
Msumprint=[]
tr=troceador()
largo_lin=80
tr.P_go(O.este,largo_lin)
Esteprint=tr.matriz
tr.P_go(O.estt,largo_lin)
Esttprint=tr.matriz
tr.P_go(list(Frec.out),largo_lin)
Frecprint=tr.matriz
tr.P_go(list(num_aa.out),largo_lin)
Numprint=tr.matriz

for i in range(len(O.Malin)):
	tr.P_go(O.Malin[i],largo_lin)
	Malinprint.append(tr.matriz)
for i in range(len(E.Msum)):
	tr.P_go(E.Msum[i],largo_lin)
	Msumprint.append(tr.matriz)
	


if "-fileout" in sys.argv:
	index=1+sys.argv.index("-fileout")
	out_name=sys.argv[index]
else:
	out_name="OUT"

salida=open(out_name,"w")
salida.write("Output of Sguallino 2.0-beta\n")
salida.write(("Line1: Index of the residue in the sequence "+text_sec+" (printed every 10 residues)\n"))
salida.write("Line 2: Relative abundance in proteins for eachr esidue in the problems sequence. \n")
salida.write("Line 3: Predicted secondary structure for the problem sequence.\n")
salida.write("Line 4: Empirical secondary structure for the first template.\n")
salida.write("Lines 5 to 5+X: alignment\n")
salida.write("Lines 6+X to 5+2X: quality of each alignment\n")
salida.write("Alignment quality key: Polarity: 4, Beta-sheet tendendcy: 2, alpha-helix tendency: 1\n")
salida.write("Abundance key: F: >6.5%, C: >5.2%, U: >3.2 %, R: =< 3.2.0%\n\n")
#salida.write("Clave de abundancia: F: >6.0%, C: >4.0%, U: >2.0%, R: =< 2.0%\n\n") #version antigua

E.identity.insert(0,100)
salida.write("Identity percentages are relative to the problem sequence.\n")
for i in range(len(O.nombres_sec)):
	salida.write(O.nombres_sec[i]+" Identity: "+str(E.identity[i])+"%\n")
salida.write("\n\n")

for i in range(len(Malinprint[0])):
	salida.write(Numprint[i])
	salida.write("\n")
	salida.write(Frecprint[i])
	salida.write("\n")
	salida.write(Esttprint[i])
	salida.write("\n")
	salida.write(Esteprint[i])
	salida.write("\n")
	for j in Malinprint:
		salida.write(j[i])
		salida.write("\n")
	for j in Msumprint:
		salida.write(j[i])
		salida.write("\n")
	salida.write("\n\n")
salida.close()
