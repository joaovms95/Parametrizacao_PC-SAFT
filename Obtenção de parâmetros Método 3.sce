//Cálculo dos parâmetros da PC-SAFT usando tabela 1 de Al Ajmi et al. 2011: F1
//Método 3
clear
clc

exec("C:\Users\jvms9\OneDrive\Área de Trabalho\Finalização de dissertação\Comparação Al-Ajmi CBTERMO versão corrigida\F1\Liang CM1\Função PC-SAFT para obter Z.sce")

//Massas molares de acordo com tabela 4.6 do Riazi (2005)), p.162
x = [4.114;4.181;3.698;3.412;2.954;2.535;2.332;1.958;1.854;1.661;1.431;1.306;1.259;1.126;0.997;0.94;0.845;0.762;0.727;0.625;0.587;0.583;0.542;0.528;0.541;0.458;0.44;0.399;0.391;6.659] //%molar do C7 a C36
x = x./100 //frações molares do C7 a C36
xn = x/sum(x) //frações molares normalizadas
MM = [95;107;121;136;149;163;176;191;207;221;237;249;261;275;289;303;317;331;345;359;373;387;400;415;429;443;457;471;485]//Massa molar C7 a C35

x7mais = sum(x)
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

M = [MM;M36mais]

//Cálculo de SG de cada MCN através da correlação de Soreide (1989)
function z = s(Cf)
    for i = 1:length(M)
        SG(i) = 0.2855 + Cf*(M(i)-66)^0.13
    end
    soma5 = 0
    for i = 1:length(M)
        soma5 = soma5 + xn(i)*M(i)/SG(i)
    end
    z = SG7mais - xn7mais*M7mais/soma5
endfunction
Cfchute = 0.5
Cf = fsolve(Cfchute,s)
SG = 0.2855 + Cf*(M-66)^0.13
//Adicionando SG e M do C6 de acordo com a tabela 4.6 do Riazi(2005) p.162
SG = [0.690;SG]
M = [84;M]

//Cálculo de propriedades como densidade, índice de refração de cada MCN
In = 0.34-exp(2.30884-2.96508*M^0.1) //Parâmetro do índice de refração p.161 - 162 Riazi (2005) 20°C
nr = sqrt((1+2*In)./(1-In)) //p. 66 Riazi(2005), índice de refração a 20°C
nr15 = nr - 0.0004*(273.15+15.5-293.15) //p. 67 Riazi (2005), índice de refração a 15,5°C
In15 = (nr15^2-1)./(nr15^2+2)
d15 = 0.999*SG //densidade a 15,5°C
d = In.*d15./In15 //densidade a 20°C p.68 Riazi(2005)

//Cálculo das frações PNA ref: p.126-127 Riazi(2005) (método ndm)
v = 2.51*(nr-1.475)-(d-0.851)
w = (d-0.851)-1.11*(nr-1.475)
for i = 1:length(SG)
    if v(i) > 0 then
        a(i) = 430
        b(i) = 0.055
    else
        a(i) = 670
        b(i) = 0.080
    end
end
S = zeros(length(v),1) //teor de enxofre
for i = 1:length(SG)
    if w(i) > 0 then
        pCr(i) = 820*w(i)-3*S(i)+10000/M(i) //percentual de carbono em aneis
    else
        pCr(i) = 1440*w(i)-3*S(i)+10600/M(i) //percentual de carbono em aneis
    end
end
pCa = a.*v + 3660./M //Percentual de carbono aromático
pCn = pCr - pCa //Percentual de carbono naftênico
pCp = 100 - pCr //Percentual de carbono parafínico
xCa = pCa/100 //Fração de carbono aromático
xCn = pCn/100 //Fração de carbono naftênico
xCp = pCp/100 //Fração de carbono parafínico

//Corrigindo frações negativas que ocorreram para SCN com M<200 (a fração negativa é descontada de todas as frações, após isso, é feita rernormalização)
PNA = [xCp,xCn,xCa]
[linh,col] = size(PNA)
for i = 1:linh
    for j = 1:col
        if PNA(i,j)<0 then
            const = PNA(i,j)
            PNA(i,1) = PNA(i,1)-const
            PNA(i,2) = PNA(i,2)-const
            PNA(i,3) = PNA(i,3)-const
            soma = sum(PNA(i,:))
            PNA(i,1) = PNA(i,1)/soma
            PNA(i,2) = PNA(i,2)/soma
            PNA(i,3) = PNA(i,3)/soma
        end
    end
end
xCp = PNA(:,1)
xCn = PNA(:,2)
xCa = PNA(:,3)

//Cálculo dos parâmetros m, epsilon_K e sigma da PC-SAFT para cada SCN (Liang et al. (2014)) CM1
for i = 1:length(SG)
    m(i) = (0.02569*M(i)+0.8709)*xCp(i) + (0.02254*M(i)+0.6827)*xCn(i) + (0.02576*M(i)+0.2588)*xCp(i)
    msigma3(i) = (1.7284*M(i)+18.787)*xCp(i) + (1.7115*M(i)+1.9393)*xCn(i) + (1.7539*M(i)-21.324)*xCp(i)
    sigma(i) = (msigma3(i)/m(i))^(1/3)
    me_K(i) = (6.8248*M(i)+141.14)*xCp(i) + (6.4962*M(i)+154.53)*xCn(i) + (6.6756*M(i)+172.40)*xCp(i)
    e_K(i) = me_K(i)/m(i)
end
