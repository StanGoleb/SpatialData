[PL]  
  
Program do przeprowadzenia prostej analizy i graficznego przestawienia danych z próbek z biopsji pacjentów rakowych.  
  
=====INSTRUKCJA=====  
1.Przed pierwszym uruchomieniem programu, dodać katalog "if_data" z danymi o próbkach z biopsji (IF1). Katalog należy stworzyć wewnątrz katalogu ze wszystkimi plikami programu.  
  
2.Głównym plikiem, ktory należy uruchomić, jest visualization.py. W tym celu użyć komendy:  

streamlit run <sciezka_do_katalogu_roboczego>/visualization.py  
  
3.Pliki clusters.py i colors.py są niezbędne do przeprowadzenia obliczeń. Zostaną automatycznie wywołane przez skrypt wizualizacji.  

  
[ENG]  
  
Simple script for analysis and graphical representation of data from biopsy samples of cancer patients  
  
=====INSTRUCTIONS=====  
  
1.Before the first run of the program, create a directory named “if_data” containing data about biopsy samples (IF1). Place this directory inside the main directory where all program files are located.  
2.The main file to execute is visualization.py. To run the script use the command:  
  
streamlit run <path_to_working_directory>/visualization.py  
  
3.The files clusters.py and colors.py are necessary for the calculations and will be automatically invoked by the visualization script.
