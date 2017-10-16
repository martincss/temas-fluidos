# temas-fluidos

La guía 2 (y próximas) contienen los datos obtenidos con ```GHOST```. Para que el solver guarde los outputs binarios directamante a esta carpeta (guia 2), en el archivo ```parameter.inp```, cambiar la linea ```odir``` con el directorio deseado.
Por ejemplo:  ```odir = "/home/usuario/Desktop/temas-fluidos/guia_2/"```

Esto solo mueve los archivos binarios (de componentes de velocidad, por ejemplo), mientras que otros outputs en formato txt se generan directamente en el directorio ```/bin```. Si nos ubicamos en ```/bin``` pueden ser fácilmente movidos con ```mv *.txt /home/usuario/Desktop/temas-fluidos/guia_2/```.

Si bien conviene guardar los archivos de datos en el directorio local del repositorio, sus formatos ```.txt``` y ```.out``` o sus carpetas deben ser agregados al ```.gitignore``` para así evitar que se suban al repo, teniendo en cuenta su tamaño intractable.
