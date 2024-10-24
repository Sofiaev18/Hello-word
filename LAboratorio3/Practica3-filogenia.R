library(ape)
library(phangorn)
library(phytools)

#Para empezar convertimos el archivo en un formato compatible con la librería phangorn, indicando el formato: “FASTA” y el tipo de secuencias: “AA” (aminoácidos)
#Asegúrate de que el archivo con las secuencias alineadas está en la carpeta donde está el Proyecto de RStudio. Y en file = “” coloca el nombre CORRECTO del archivo de las secuencias Ejemplo: “fraxatin_aligned.fasta”
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                        format = "FASTA", type = "AA")
fraxatin

#Ahora lo que se hace es crear una matriz de distancia para poder crear árboles de distancia o parsimonia a través de ella. En este caso la clase debe ser AAbin (Amino Acid Sequences), por eso transformamos el objeto fraxatin en este tipo de clase. La función dist.aa (método o función de objeto o clase) calcula una matriz de distancias por pares de las secuencias de aminoácidos partir de una objecto de clase AAbin utilizando un modelo de evolución de aminoácidos (por ejemplo Dayhoff).
#Estos árboles son árboles de similaridad, no representan inferencia evolutiva.
matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist

#Con la matriz de distancia perdemos los caracteres a favor de las diferencias en caracteres entre especies. Cuando aparece un valor de 0, significa que no hay diferencia al nivel de los caracteres (aminoácidos) entre las sequencias de dos especies. Ahora creamos un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) usando la matriz de distancia que acabamos de calcular.
#Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)

#Ahora hacemos un árbol con el método de unión de vecinos (NJ) usando la misma matriz de distancias:
arbolNJ <- nj(matrizdist)
plot(arbolNJ)

#Este último árbol puede ser distinto del árbol creado con el método UPGMA. Para personalizar los árboles podemos agregar argumentos a parámetros como cex, para el tamaño de la letra, edge.color, para el grosos de las ramas, etc. También se puede escoger entre diferentes visualizaciones de árbol como filograma, cladograma, radial y demás.
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)

plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)

#Además de plot podemos graficar árboles con el método plotTree del paquete phytools, el cual es compatible con ape y con phangorn
plotTree(arbolNJ)

#A este también le podemos hacer modificaciones.
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)


#En los árboles, sin cambiar la topología, se puede cambiar el orden en que los grupos se ven Por ejemplo, se pueden ordenar las puntas de manera alfabética (en la medida de lo posible), o con los grupos más derivados hacia uno de los lados del árbol. Para escalerizar hacia la derecha:
plotTree(ladderize(arbolNJ))

#Para guardar un árbol se usa el comando: 
write.tree(arbolNJ, file = "file_name.nex").

#se lee con:
read.tree(file = "file_NJ.nex").


#ENRAIZAR
#Hasta ahora los árboles no estaban enraizados, para ello se usa la función root del paquete ape.para poner por raíz las secuencias de fraxatina de Ornitorinco pasamos por argumento el nombre de la secuencia que corresponde a esta secuencia en el parámetro outgroup (Ornitorrinco).
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)

#Se puede hacer lo mismo con el árbol de UPGMA
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)

#Para ver ambos a la vez
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)


#ÁRBOLES DE PARSIMONIA

#La parsimonia busca disminuir el número de pasos que explican un árbol evolutivo contando el número de cambios de cada uno de los caracteres
#En los métodos de distancia UPGMA y NJ se llega a un árbol pero con parsimonia se evalúan varios árboles 
#La parsimonia no usa todos los aminoácidos de la secuencia. Se eliminan los que son constantes en todos los taxones o que son variables pero no proporcionan información, Un caracter es informativo si tiene al menos dos estados de caracter y por lo menos dos de estos estados ocurren con una frecuencia mínima de dos.
#Para estimar árboles de máxima parsimonia existen varias posbilidades, la más sencilla es partir de árboles de distancia. Se utiliza un árbol de inicio obtenido por distancia y se cuenta su número de pasos
  #EJEMPLO--> número de pasos del árbol arbolUPGMAraiz.
parsimony(arbolUPGMAraiz, fraxatin) #313 --> este árbol tiene 313 pasos, aunque tenga raíz o no debe tener los mismos pasos 
  #EJEMPLO --> número de pasos del árbol pero sin raiz 
parsimony(arbolUPGMA, fraxatin)  #313 


#Con el método optim.parsimony se obtiene el árbol con mejor parsimonia. Este método permite encontrar árboles bajo máxima parsimonia usando árboles de distancia de inicio.
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin). #Puntuación final es 307 despues de 2nni operaciones 
  #Con el árbol NJ
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin). #puntuación final es 307 después de 1nni operaciones 

#Hacer parsimonia con un número de taxones muy alto se convierte en una labor computacional ardua. Más allá de 10 taxones hay un número de combinaciones difíciles de obtener computacionalmente y se deben hacer búsquedas heurísticas, esto significa que no se evalúan todas las posibilidades de árboles.
#Los métodos heurísticos buscan algunos árboles seleccionando un árbol inicial y cortando ramas y poniéndolas en otro lugar del árbol (branch swapping). Se puede interpretar como una muestra representativa de una problación para hacer inferencias a partir de ella. Hay varias formas de búsquedas heurísticas como:

#Otra estrategia para hacer el proceso de búsqueda de árbol con mayor parsimonia es con el algoritmo de búsqueda pratchet. El cual tiene los siguientes pasos:
    #1. Generar un árbol de inicio con algún nivel de intercambio de ramas 
    #2. Seleccionar al azar un subconjunto de caracteres (aminoácidos) y darles más peso, esto quiere decir que se usan con mayor frecuencia dichos caracteres para hacer los análisis. La cantidad de caracteres que se seleccionan es establecida por el usuario, típicamente entre 5 y 25% de los caracteres.
    #3. Realizar intercambio de ramas (branch swapping).
    #4. Iterar pasos 1 a 3 entre 50 y 200 veces. El algoritmo anterior es computacionalmente más amigable y rápido. Probémoslo usando la función pratchet:
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
fraxatin_parsimonia. #4 phylogenetic trees 


#Para poderlos comparar es necesario enraizarlos.
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)

#Para escoger el sólo un árbol con igual parsimonia hay que conseguir un árbol de consenso.
#Para hacer un árbol de consenso estricto podemos usar el método ape con parámetro p de 1, que corresponde a un 100% de consenso entre ramas.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)

#Para un árbol menos estricto podemos cambiar el valor del parámetro p:
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)


#BOOTSTRAP
#Para dar soporte a los árboles se puede hacer una serie de seudoréplicas con remplazamiento (bootstrapping) de una matriz
#Con cada réplica se hace un árbol consenso; la veces que un grupo se repita en el conjunto de réplica es el valor de soporte del nodo. Esta técnica, ampliamente utilizada, proporciona evaluaciones de «confianza» para cada clado de un árbol observado, basándose en la proporción de árboles bootstrap que muestran ese mismo clado. Normalmente un nodo con un soporte mayor a 80% es un buen soporte.
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)

#En este caso usamos la función pratchet, si usáramos otra se demoraría más. También usamos un número de réplicas irrisosio pues este es un ejercicio demostrativo y no es necesario gastar mucho tiempo en él. La rutina anterior genera entonces 10 árboles pseudoréplicas.
plot(arbolesbootstrap, cex = .6)

#Ahora se genera un consenso del 60%
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

#MODELOS PROBABILÍSTICOS 
#Los caracteres en este caso no son sólo caracteres si no que también representan dinámicas evolutivas. A diferencia de la parsimonia, con estos modelos todos los caracteres son útiles, incluso los que son constantes entre todos los taxones

#ÁRBOLES DE MÁXIMA VEROSIMILITUD 
#Para calular la verosimilitud de un árbol segun un alineamiento de secuencias se usa un modelo de sustitución de aminoácidos. También se tienen en cuenta la frecuencia de de ocurrencia de los aminoácidos
#Para ello hay una serie de pasos: 
    #1. Se genera un árbol enraizado de inicio de cualquier tipo, puede ser al azar.
    #2. Se calcula la verosimilitud de casa sitio (de aminoácidos) usando un árbol.
    #3. Se calcula la verosimilitud total del árbol por nodo, se consideran todos los escenario posibles que pudieron dar origen a los estados de caracter observados para dicho sitio. La verosimilitud de un sitio es, entonces, la suma de verosimilitudes de la reconstrución de este sitio (en forma de árbol) dado un modelo de evolución o sustitución.

#Para el ejemplo usamos de nuevo nuestro objeto fraxatin de clase phy y creamos un árbol al azar de 11 ramas (porque tenemos 11 secuencias) con rtree como punto de partida.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)

#En seguida lo enraizamos por las secuencias de Ornitorinco para poderlo visualizar mejor. Además los «escalerizamos» hacia la derecha y le agregamos escala; aquí la longitud de la rama sí es significativa, indica cantidad de cambio en cuanto a sustituciones de aminoácidos.
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()

#A partir del arbol de arriba se puede iniciar la búsqueda del mejor árbol por máxima verosimilitud. Lo primero que se hace es carcular la verosimilitud del árbol dadas las secuencias. Con pml (Phylogenetic maximum likelihood), podemos computar tal verosimilitud.
ajustado <- pml(arbolazarR, fraxatin).      
ajustado

    #REsultado ### model: Mk 
## loglikelihood: -4367.895 
## unconstrained loglikelihood: -1479.871 
## Rate matrix:


#La información que tiene el objeto ajustado nos reporta la verosimilitud del árbol al azar que habíamos creado, que es -4348.064. También reporta un modelo de substitución general, el cual tal vez no se ajuste bien a los datos. Lo que hay que hacer es encontrar un árbol que optimice la verosimilitud usando un modelo de sustitución; para esto vamos a usar el método optim.pml del paquete phangorn, el cual computa la verosimilitud de un árbol filogenético dado un alineamiento múltiple de secuencias y un modelo de evolución de AA. Toma como argumentos un objeto de clase pml, el tipo de modelo que se quiere usar así como el tiempo de rearreglo para los árboles.
ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

#Para ver el árbol oculto usamos $tree. También lo enraizamos.
ajustadoconDay$tree

ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()

#El árbol anterior fue generado usando la matriz de sustitución de Dayhoff. Pero se pueden usar diferentes modelos.
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")

#Podemos comparar los modelos calculando el Criterio de información de Akaike AIC:
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)

#La primera columna corresponde a los grados libertad. Según el criterio anterior, el mejor modelo que se ajusta con los datos (con el AIC 
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")

mejorarbol

mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()

#ESTE ES EL MEJOR ÁRBOL SIGUIENDO EL MODELO DE MÁXIMA VEROSIMILITUD 


















