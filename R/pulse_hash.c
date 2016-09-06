//
// pulse_hash.c
//
// 
//
//
//
//
//

typedef struct pulse_container {

	double *time[];
	double theta[][2];
	double mean_contrib[]fitlength[];
	double death_rate;
	double total_death_rate;

} pulses; 


string *deleteEntry(string *dynamicArray, int &size, string entryToDelete)
{// create a new dynamic array 1 element larger than dynamicArray
   string *newArray = new string[size - 1];
 
   // copy all elements from dynamicArray into new array
   for(int i = 0; i < size; i++)
   {
       if(dynamicArray[i]!=entryToDelete)
           newArray[i] = dynamicArray[i];
   }
 
   size--;
 
   // delete dynamicArray
   delete [] dynamicArray;
 
   // and return the new array
   return newArray;
}
 
    
 
string *addEntry(string *dynamicArray, int &size, string newEntry)
{
   // create a new dynamic array 1 element larger than dynamicArray
   string *newArray = new string[size + 1];
 
   // copy all elements from dynamicArray into new array
   for(int i = 0; i < size; i++)
   {
       newArray[i] = dynamicArray[i];
   }
 
   // add the new entry onto the end of the new array
   newArray[size] = newEntry;
 
   // increment size
   size++;
 
   // delete dynamicArray
   delete [] dynamicArray;
 
   // and return the new array
   return newArray;
}
