//Cálculo dos parâmetros da PC-SAFT usando tabela 1 de Al Ajmi et al. (2011): F1
//Método 1
clear
clc

exec("C:\Users\jvms9\OneDrive\Área de Trabalho\Comparação Al-Ajmi CBTERMO\F1\14nc-4, hosseinifar, assareh\Função PC-SAFT para obter Z.sce")

//Massas molares
x = [4.114;4.181;3.698;3.412;2.954;2.535;2.332;1.958;1.854;1.661;1.431;1.306;1.259;1.126;0.997;0.94;0.845;0.762;0.727;0.625;0.587;0.583;0.542;0.528;0.541;0.458;0.44;0.399;0.391;6.659] //%molar do C7 a C36
x = x./100 //frações molares do C7 a C36
xn = x/sum(x) //frações molares normalizadas
NC = []
for i = 1:29
    NC = [NC;[6+i]]
end
MM = 14*NC-4//Massa molar C7 a C35

M7mais = 291
SG7mais = 0.8945/0.999
xn7mais = sum(xn)//fração molar 7+ normalizada
xn36mais = xn(length(xn)) //fração molar 36+ normalizada

//Cálculo de MM36+
soma = 0
for i = 1:length(xn)-1
    soma = soma+MM(i)*xn(i)
end
M36mais = (M7mais-soma)/xn36mais

//Concatenando os vetores para obter todos os componentes
M = [14*6-4;MM;M36mais] //o primeiro termo é do C6 (Riazi (2005), tabela 4.2)

//Cálculo de SG de cada SCN através da correlação de Hosseinifar et al. (2014)
A = 0.0089
B = 0.70799
C = -0.14064
D = 0.13199
SG = (A*M.^B+C).^D

//Cálculo dos parâmetros m, epsilon_K e sigma da PC-SAFT para cada SCN (Assareh et al. (2016))
for i = 1:length(M)
    m(i) = 33.58+0.08816*M(i)-90.75*SG(i)-0.07727*M(i)*SG(i)+61.01*SG(i)^2
    me_K(i) = 3372+11.24*M(i)-8955*SG(i)-5.925*M(i)*SG(i)+6136*SG(i)^2
    e_K(i) = me_K(i)/m(i)
    msigmacubo(i) = -75.14+2.848*M(i)+231.7*SG(i)-1.288*M(i)*SG(i)-186.9*SG(i)^2
    sigma(i) = (msigmacubo(i)/m(i))^(1/3)
end
