*********************************************************************
*                      Ejemplo de Estimacion usando                 *
*                             Modulo sp de stata                    *
*              Generando y exportando diferentes  matrices          *
*                             Noviembre 2021                        *
*********************************************************************

set more off
clear all
clear mata

global rutabox "C:\Users\lmoralzu\Dropbox\"
*global rutabox "C:\Users\zurit\Dropbox\"

*Open shape in Stata
/*cd "/Users/Eleo 1/Dropbox/Papers/Papers CI/Spatial_Analisis_IC/Estimation 3paperv3"*/
*cd "E:\Spatial_Analisis_IC\Estimation 3paperv3"
*cd "C:\Users\zurit\Dropbox\Spatial_Analisis_IC\LM\Estimation 3paperv3"
 cd "${rutabox}\Spatial_Analisis_IC\LM\Estimation 3paperv4"   /*06/05/2022->Cambió una variable, copio y pego la nueva carpeta en LM*/

 
spshape2dta mpio, replace
*use  mpio.dta
*use  mpio_shp.dta
 

*Database
use "mpio.dta", clear
spset

list _ID _CX _CY if NOMBRE_MPI=="MEDELLIN" | NOMBRE_MPI=="CALI"
/*Distancias estan en metros*/
grmap, activate
grmap HECTARES

/*Creando Matrices con el comando spmatrix create*/
use "mpio.dta", clear
spset
spmatrix create contiguity W, replace

/*Pareando con el Panel de la Estimacion*/
use "${rutabox}\Spatial_Analisis_IC\LM\Stata_CI_Paper3\dataset for paper3_v4.dta", clear  /*06/05/2022*/
rename id idyear
duplicates report idmun year
gen MPIOS = idmun

gen idmunN=idmun
destring idmunN, replace

xtset idmunN year
spbalance
merge m:1 MPIOS using "mpio.dta", keep(3) nogen
spset

keep if year==2001

* * * Creando y exportando diferentes matrices * * *
*Creando Matrices Vecindad*
spmatrix create contiguity Wq1,  first              normalize(row) replace
spmatrix create contiguity Wq2,  first second       normalize(row) replace
spmatrix create contiguity Wq3,  first second(0.5) normalize(row) replace

*Creando Matrices Distancia*
forval d=1/5{
local r=1/(10*1000*`d')
dis(`r')
spmatrix create idistance  W`d'k, vtruncate(`r') normalize(row) 
}

*Matriz de distancia normalita*
spmatrix create idistance  WDk,  normalize(row)

spmatrix dir

*Exportando*
global Wmatrix "Wq1 Wq2 Wq3 WDk W1k W2k W3k W4k W5k"
foreach W of global Wmatrix{
spmatrix export `W'  using "${rutabox}\Spatial_Analisis_IC\LM\Stata_CI_Paper3\\`W'.txt" , replace
spmatrix export `W'  using "${rutabox}\Spatial_Analisis_IC\LM\Stata_CI_Paper3\\`W'.csv" , replace
}















