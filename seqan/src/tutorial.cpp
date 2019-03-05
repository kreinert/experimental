#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/find.h>



using namespace seqan;
using namespace std;

int main()
{
    // One possible solution to the first sub assignment
    String<AminoAcid> text = "VQLPPPGQQLQVLLPKKPP";
    Index<String<AminoAcid>, FMIndex<> > index(text);
    
    StringSet<String<AminoAcid> > stringSet;
    typedef Iterator<StringSet<String<AminoAcid> > >::Type TStringSetIterator;
    
    appendValue(stringSet, "VQLPPK");
    appendValue(stringSet, "KKP");
    appendValue(stringSet, "QVLL");
    
    Finder<Index<String<AminoAcid>, FMIndex<> > > myFinder(index);
    
    TStringSetIterator it;
    for (it = begin(stringSet); it != end(stringSet); ++it)
    {
        std::cout << *it << ":\n";
        while (find(myFinder,*it)) {
            cout << "Found " << *it << " at " << position(myFinder) << endl;
        }

        clear(myFinder);
    }

 
  
    return 0;
}
