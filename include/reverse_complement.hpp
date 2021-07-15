#ifndef REV_COMP_HPP
#define REV_COMP_HPP

string reverse_complement(string sequence){
  reverse(sequence.begin(), sequence.end());
  string reverseComp;
  for(int i=0; i < sequence.size(); ++i){
    char letter = sequence[i];
    if(letter == 'A'){
      reverseComp.append("T");
    }
    if(letter == 'T'){
      reverseComp.append("A");
    }
    if(letter == 'C'){
      reverseComp.append("G");
    }
    if(letter == 'G'){
      reverseComp.append("C");
    }
    if(letter == 'N'){
      reverseComp.append("N");
    }
  }
  return reverseComp;
}  

#endif
