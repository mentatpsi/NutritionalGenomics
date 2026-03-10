# NutritionalGenomics
This Flask Package is part of the OSGenome Ecosystem and features a way to learn more about your individual nutritional needs

It leverages my OSGenome project's processing on 23andMe raw data and produces a table on SNPs associated with results on the Ketogenic diet as well as other nutritional guidelines. Including impact on BMI with Saturated Fat, Omega 3s, impact on Folate Metabolism, and more. It handles flipping genotypes based on orientation. The SNPs associated with the ketogenic diet were previously researched and greatly expanded on to produce 18 specific nutritional genotypes on my own raw data. I was previously doing this on Jupyter. It features Python's Flask and Pandas for the back end and some basic CSS and HTML, all crafted by AI. Also provides a link to the SNPedia page of the RSID for more details.

To learn more:
https://github.com/mentatpsi/OSGenome

Once you've set up the snpDict.json with OSGenome, place it in the directory you cloned the project in to learn more about how genetics can help guide nutrition. Don't forget to run pip install -r requirements.txt to get required dependencies.


## Example
![Example of Kendo Grid](https://github.com/mentatpsi/NutritionalGenomics/blob/main/screenshots/NutritionalGenomics1.png)
