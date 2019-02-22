#include <iostream>
#include <vector>
#include <algorithm>

int main(){

    double sumEnergy, energy, J;
    double Emin = 9999., Emax = -9999.;
    std::vector<int> lattice;
    std::vector<std::vector<int> > storedLattice;
    int row, col, start, end;
    std::cout << "Rows: ";
    std::cin >> row;
    std::cout << "Cols: ";
    std::cin >> col;
    std::cout << "J: ";
    std::cin >> J;
 //   std::cout << "start";
 //   std::cin >> start;
    std::cout << "end";
    std::cin >> end;  
    //creates the first empty lattice
    for (int i = 0 ; i < row * col; i++){
        lattice.push_back(0);
    }
    storedLattice.push_back(lattice);

    //creates all the possible combinations in one big chunk
    for (int j = end - 1; j >= 0; j--){
        lattice[j] = 1;
	if (end == (row*col)){
		storedLattice.push_back(lattice);
	   	while(std::next_permutation(lattice.begin(), lattice.end())){
     			 storedLattice.push_back(lattice);
        	}
	}
    }
    if( end != (row*col)){
	    for(int j = end - 1; j>= 0; j--){
		storedLattice.push_back(lattice);
		while(std::next_permutation(lattice.begin(), lattice.end())){
			 storedLattice.push_back(lattice);
		}
	    }
    }

    for(int i = 0; i < storedLattice.size(); i++){
        int counter = 0;
        sumEnergy = 0;
        //for energy in rows
        do {
            sumEnergy += (storedLattice[i][counter * col] * storedLattice[i][counter * col + (col-1)]);
            for (int j = 0; j < col - 1; j++) {
                int k1, k2;
                k1 = counter * col + j;
                k2 = k1 + 1;
                energy = storedLattice[i][k1] * storedLattice[i][k2];
                sumEnergy += energy;
            }
            counter++;
        }while ( counter < row);

        counter = 0;
        //for energy in cols
        do {
            sumEnergy += (storedLattice[i][counter * row] * storedLattice[i][counter * row + (row-1)]);
            for (int j = 0; j < row - 1; j++) {
                int k1, k2;
                k1 = counter * col + j;
                k2 = k1 + 1;
                energy = storedLattice[i][k1] * storedLattice[i][k2];
                sumEnergy += energy;
            }
            counter++;
        }while ( counter < col);

        sumEnergy = J * sumEnergy;
        if (sumEnergy > Emax) {
            Emax=sumEnergy;
        }
        if (sumEnergy < Emin) {
            Emin=sumEnergy;
        }
    }

    std::cout << "total number of states is: " << storedLattice.size() << std::endl;
    std::cout << "Emax is: " << Emax << std::endl;
    std::cout << "Emin is: " << Emin << std::endl;


    //prints all the possible combinations
    for (int l = 0; l < storedLattice.size(); l++){
        for (int i = 0 ; i < row; i++){
            for(int j = 0; j < col; j++){
                int k = i * col + j;
                std::cout << storedLattice[l][k];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}