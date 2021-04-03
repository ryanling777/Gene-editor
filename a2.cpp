/*Gene Editing Program by Ryan Ling 25/3/2021.
This code is able to receive FASTA (.fna) files and is able to conduct the functions of finding, adding, deleting and replacing sequences
in multiple files by using custom linked list structures and Object Orientated Programming concepts. */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "tictoc.h"

using namespace std;

void check_whitespace(vector<string> &filenames); //function to remove leading and ending white spaces.

//each Node stores a character.
class Node {
private:
  Node* p_next;
  Node* p_previous;
  char base_pair;
public:
  Node();
  Node(char basepair, Node* p_next_in, Node* p_previous_in);
  char get_basepair();
  Node* get_next();
  Node* get_previous();
  void set_next(Node* p_next_in);
  void set_previous(Node* p_previous_in);
  void set_basepair(char basepair);
};

Node::Node() {
  p_next = nullptr;
  p_previous = nullptr;
  base_pair = 'A';
}

Node::Node(char basepair, Node* p_next_in, Node* p_previous_in) {
  p_next = p_next_in;
  p_previous = p_previous_in;
  base_pair = basepair;
}

char Node::get_basepair() {
  return base_pair;
}

Node* Node::get_next() {
  return p_next;
}

Node* Node::get_previous() {
  return p_previous;
}

void Node::set_next(Node* p_next_in) {
  p_next = p_next_in;
}

void Node::set_previous(Node* p_previous_in) {
  p_previous = p_previous_in;
}

void Node::set_basepair(char basepair) {
  base_pair = basepair;
}




//Stores the sequences from each file from the user. Each Doublylinkedlist = 1 file sequence.
class Doublylinkedlist {
private:
  Node* p_head;
  Node* p_tail;
  Doublylinkedlist* next_ll; //points to the next linked list
  Doublylinkedlist* previous_ll; //points to previous linked list.
public:
  Doublylinkedlist();
  void append(char basepair);
  void append(Node* nodetobeadded); //ovreloaded append.
  void insert(Doublylinkedlist* original_list, Doublylinkedlist* added_sequence, int position);
  void delete_seq(int bpposition, int bplength, Doublylinkedlist* original_list);
  void replace_seq(int bpposition, int bplength, Doublylinkedlist* original_list, Doublylinkedlist* replacingseq);
  void search(Doublylinkedlist* desired_sequence, Doublylinkedlist* original_list);
  void printlist();
  Doublylinkedlist* get_nextll();
  void set_nextll(Doublylinkedlist* p_next_in);
  void pop_back();
  Node* get_p_head();
  Node* get_p_tail();
  Doublylinkedlist* get_previousll();
  void set_previousll(Doublylinkedlist* p_previous_in);
  int size(); //returns size of sequence (linked list)
  void set_p_head(Node* p_head_in);
  void set_p_tail(Node* p_tail_in);
  ~Doublylinkedlist();

};

Doublylinkedlist::Doublylinkedlist() {
  p_head = nullptr;
  p_tail = nullptr;
  next_ll = nullptr;
  previous_ll = nullptr;
}

void Doublylinkedlist::append(char basepair) {

  if (p_head == nullptr) {
    //if linked list is empty
    p_head = new Node(basepair, p_head, p_head);
    p_tail = p_head;

  }

  else {

    Node* p_seek = p_tail;
    Node* new_node = new Node(basepair, 0, p_seek); //sets p next to null and p previous to pseek
    p_seek->set_next(new_node);
    p_tail = new_node;

    }

}

void Doublylinkedlist::append(Node* nodetobeadded) {
  if (p_head == nullptr) {
    //if linked list is empty.
    p_head = nodetobeadded;
    p_tail = p_head;


  }
  else {

    nodetobeadded->set_previous(p_tail);
    p_tail->set_next(nodetobeadded);
    p_tail = nodetobeadded;

    }

}

void Doublylinkedlist::printlist() {
  for (Node* p_seek = p_head; p_seek != nullptr; p_seek = p_seek->get_next()) {
    cout << (p_seek->get_basepair());
  }
}

Doublylinkedlist* Doublylinkedlist::get_nextll() {
  return next_ll;
}

void Doublylinkedlist::set_nextll(Doublylinkedlist* p_next_in) {
  next_ll = p_next_in;
}

void Doublylinkedlist::pop_back() {
  Node* p_temp;

  //if there is nothing insidee the list
  if (p_tail == nullptr) {
    cout << "there is nothing currently in the list to pop" << endl;
  }
  //if there is only one node remaining.
  else if (p_head == p_tail) {
    p_temp = p_head;
    p_head = nullptr;
    p_tail = nullptr;
    delete p_temp;
  }
  else {
    p_temp = p_tail;
    p_tail = p_tail->get_previous();
    p_tail->set_next(nullptr);
    delete p_temp;
  }


}

Node* Doublylinkedlist::get_p_head() {
  return p_head;
}

Node* Doublylinkedlist::get_p_tail() {
  return p_tail;
}

Doublylinkedlist* Doublylinkedlist::get_previousll()  {
  return previous_ll;
}
void Doublylinkedlist::set_previousll(Doublylinkedlist* p_previous_in) {
  previous_ll = p_previous_in;
}


void Doublylinkedlist::insert(Doublylinkedlist* original_list, Doublylinkedlist* added_sequence, int position) {

  Node* p_seek = original_list->get_p_head();
  Node* p_temp; // both are temporary pointers to store certain positions without changing where p_seek is pointing to.
  Node* p_temp2;
  int added_seqlength = added_sequence->size(), num_previous = 0;
  char previous[10]; //stores previous 10 characters.

  cout << "Base pair positions: [" << position << ":" << position + added_seqlength-1 << "]" << endl;
  cout << "Base pair length: " << added_seqlength << endl;
  //if user wants to add to the start.
  if (position == 0) {

    cout << "Prev 10 pairs: " << 0 << endl;
    cout << "Region of interest: ";
    added_sequence->printlist();

    //connects end of user sequence to head of original sequence.
    (added_sequence->get_p_tail())->set_next(original_list->get_p_head());
    p_seek->set_previous(added_sequence->get_p_tail());

    cout << endl;
    cout << "Next 10 pairs: ";

    p_temp = original_list->get_p_head(); //for printing next 10 values.
    for (int i = 0; i < 10; i++) {
      if (p_temp != nullptr) {
        cout << p_temp->get_basepair();
        p_temp = p_temp->get_next();
      }
    }
    cout << endl;

    original_list->set_p_head(added_sequence->get_p_head()); //assing p head to front of user sequence (that is now connected to main sequence)

  }

  //if user wants to add at the end of the original list
  else if (position == original_list->size()) {

    p_temp = original_list->get_p_tail();
    cout << "Prev 10 pairs: ";
    //prints previous 10 pairs.
    for (int i = 9; i >= 0; i--) {
      if (p_temp != nullptr) {
        previous[i] = p_temp->get_basepair();
        p_temp = p_temp->get_previous();
        num_previous++;
      }

    }
    for (int i = 10-num_previous; i < 10; i++) {
      cout << previous[i];
    }

    cout << endl;
    cout << "Region of interest: ";
    added_sequence->printlist();

    //connects the end of original sequence to the head of the user sequence, and assigns tail of original sequence to end of user sequence.
    (original_list->get_p_tail())->set_next(added_sequence->get_p_head());
    (added_sequence->get_p_head())->set_previous(original_list->get_p_tail());
    original_list->set_p_tail(added_sequence->get_p_tail());

    cout << endl;
    cout << "Next 10 pairs: " << 0 << endl;

  }
  else {

    cout << "Prev 10 pairs: ";

    //p_seek traverses until the desired position.
    for (int i = 0; i < position-1; i++) {
      if (position - 1 - i < 10){
        cout << p_seek->get_basepair();
      }
      p_seek = p_seek->get_next();
    }
    cout << p_seek->get_basepair(); //for the last node just before the index to be added.
    cout << endl;
    cout << "Region of interest: ";
    added_sequence->printlist();

    p_temp = p_seek->get_next(); //p_temp points to first node after the new sequence is to be added.
    p_seek->set_next(added_sequence->get_p_head()); //the node just before insertion now points to the head of the added sequence.
    (added_sequence->get_p_head())->set_previous(p_seek); //p_previous of the head node of the added sequence points to the node just before the insertion in the original list.
    (added_sequence->get_p_tail())->set_next(p_temp); //makes p_next of the last node of the added_sequence point at the p_temp.
    p_temp->set_previous(added_sequence->get_p_tail()); //makes p_previous of first node after inserted list point back to tail of added sequence.

    cout << endl;
    cout << "Next 10 pairs: ";
    for (int i = 0; i < 10; i++) {
      if (p_temp != nullptr) {
        cout << p_temp->get_basepair();
        p_temp = p_temp->get_next();
      }
    }
    cout << endl;

  }

}

//returns size of linked list
int Doublylinkedlist::size() {
  int count = 0;
  for (Node* p_seek = p_head; p_seek != nullptr; p_seek = p_seek->get_next()) {
    count++;
  }
  return count;
}

void Doublylinkedlist::set_p_head(Node* p_head_in) {
  p_head = p_head_in;
}

void Doublylinkedlist::set_p_tail(Node* p_tail_in) {
  p_tail = p_tail_in;
}

//deletes user specified sequence.
void Doublylinkedlist::delete_seq(int bpposition, int bplength, Doublylinkedlist* original_list) {
  Node* p_seek = original_list->get_p_head(); //goes up to the position before which the deletion occurs
  Node* p_seek2 = original_list->get_p_head(); //goes up to the node just after which the deletion occurs.
  int flag = 0; // flag becomes 1 when p_seek2 passes p_seek to allow us to print out the deleted sequence.
  int original_size = original_list->size();
  //for the case where the user requests to delete a range that is more than length of sequence.
  if (bpposition + bplength > original_size) {
    cout << "Base pair positions: [" << bpposition << ":" << original_size - 1 << "]" << endl;
    cout << "Base pair length: " << -bpposition + original_size << endl;
  }

  else {
    cout << "Base pair positions: [" << bpposition << ":" << bpposition + bplength-1 << "]" << endl;
    cout << "Base pair length: " << bplength << endl;
  }

  //if user wants to add to the start.
  if (bpposition == 0) {

      cout << "Prev 10 base pairs: " << 0 << endl;
      cout << "Region of interest: ";
      for (int i = 0; i < bplength; i++) {
        if (p_seek != nullptr) {
          cout << p_seek->get_basepair();
          p_seek = p_seek->get_next(); //p_seek ends at the first head node after deletion unless it has reached the end of the sequence.
        }
      }

      cout << endl;
      //for the unique situation where user deletes everything.
      if (p_seek == nullptr) {
        original_list->set_p_head(nullptr);
        original_list->set_p_tail(nullptr);
        cout << "Next 10 base pairs: " << 0 << endl;
      }

      else {

        original_list->set_p_head(p_seek);
        p_seek -> set_previous(0);
        cout << "Next 10 base pairs: ";
        for (int i = 0; i < 10; i++) {
          if (p_seek != nullptr) {
            cout << p_seek->get_basepair();
            p_seek = p_seek->get_next();
          }
       }
        cout << endl;
      }

  }

  //for if the user inputs a position and length that deletes more than the length of the sequence.
  else if ((bpposition + bplength) >= original_list->size()) {

    cout << "Prev 10 base pairs: ";
    for (int i = 0; i < bpposition-1; i++) {
      if ((bpposition - i) <= 10) {
        cout << p_seek->get_basepair();
      }
      p_seek = p_seek -> get_next(); //p_seek finishes at node before first node to be deleted.

    }
    cout << p_seek->get_basepair(); //for the final letter before the node to be deleted.
    cout << endl;
    cout << "Region of interest: ";
    for (p_seek2 = p_seek->get_next(); p_seek2 != nullptr; p_seek2=p_seek2->get_next()) {
      cout << p_seek2->get_basepair();
    }

    p_seek->set_next(0); //removes link to the end sequence and makes p_seek the tail.
    original_list->set_p_tail(p_seek);
    cout << endl;
    cout << "Next 10 base pairs: " << 0 << endl;

  }


  else {

    cout << "Prev 10 base pairs: ";
    for (int i = 0; i < bpposition-1; i++) {
      if ((bpposition - i) <= 10) {
        cout << p_seek->get_basepair();
      }
      p_seek = p_seek -> get_next(); //p_seek finishes at node before first node to be deleted.

    }
    cout << p_seek->get_basepair(); //for the final letter before the node to be deleted.
    cout << endl;
    cout << "Region of interest: ";
    for (int i = 0; i < bpposition + bplength; i++) {
      if (p_seek2 != nullptr) {
        p_seek2 = p_seek2 -> get_next();
        if (flag == 1 && (bpposition + bplength) - i > 1) { //prints all that are deleted. Does not print the last loop as that is the first node after the deleted sequence.
          cout << p_seek2->get_basepair();
        }
        if (p_seek2 == p_seek) { //when p_seek is equal to p_seek2, flag is raised to indicate that printing can start.
          flag = 1;
        }
      }
    }


      //p_seek2 ends at the first node after deletion. No p_seek and p_seek2 are connected.
      p_seek->set_next(p_seek2);
      p_seek2->set_previous(p_seek);

      cout << endl;
      cout << "Next 10 base pairs: ";
      for (int i = 0; i < 10; i++) {
        if (p_seek2 != nullptr) {
          cout << p_seek2->get_basepair();
          p_seek2 = p_seek2->get_next();
        }

      }

    cout << endl;
  }

  cout << endl;
  cout << "The DNA sequence has been deleted" << endl;
  cout << endl;
  cout << "DNA sequence deletion result:" << endl;
  cout << "Base pair position: " << bpposition << endl;
  cout << "Previous 10 pairs: ";
  p_seek = original_list->get_p_head();

  if (bpposition == 0){ //if start of the list was deleted.
    cout << 0 << endl;
  }
  else {
    //prints previous 10 pairs.
    for (int i = 0; i < bpposition-1; i++) {
      if ((bpposition - i) <= 10) {
        cout << p_seek->get_basepair();
      }
      p_seek = p_seek -> get_next();
    }
    cout << p_seek->get_basepair();
    p_seek = p_seek->get_next();
    cout << endl;
  }
  cout << "Next 10 pairs: ";
  //if the deleted part of the sequence was at the end of the original sequence.
  if (p_seek == original_list->get_p_tail() || p_seek == nullptr) {
    cout << 0 << endl;
  }

  else {
    for (int i = 0; i < 10; i++) {
      if (p_seek!= nullptr) {
        cout << p_seek->get_basepair();
        p_seek = p_seek->get_next();
      }
    }
  }

  cout << endl;

}

void Doublylinkedlist::replace_seq(int bpposition, int bplength, Doublylinkedlist* original_list, Doublylinkedlist* replacingseq) {
  char previous[10]; //stores up to previous 10 values.
  int num_previous = 0; //counts number of previous elements.
  int replacing_length = replacingseq->size();
  Node* p_seek = original_list->get_p_head(); //goes up to the position before the replacement occurs.
  Node* p_seek2 = original_list->get_p_head(); //goes up to the node just after the replacement
  Node* p_temp; //for printing preceding values without changing p_seek.

  cout << "Base pair positions: [" << bpposition << ":" << bpposition + replacing_length-1 << "]" << endl;
  cout << "Base pair length: " << replacing_length << endl;

  //if user wants to replace the start of the sequence.
  if (bpposition == 0) {
    for (int i = 0; i < bplength; i++) {
      if (p_seek != nullptr) {
        p_seek = p_seek->get_next(); //p_seek ends at the first head node after deletion.
      }
    }

    cout << "Prev 10 base pairs: " << 0 << endl;
    cout << "Region of interest: ";

    replacingseq->printlist();
    cout << endl;
    //for the unique situation where user wants to replace the entire sequence.
    if (p_seek == nullptr) {
      original_list->set_p_head(replacingseq->get_p_head());
      original_list->set_p_tail(replacingseq->get_p_tail());
      cout << "Next 10 base pairs: " << 0 << endl;
    }

    else {
      p_seek->set_previous(replacingseq->get_p_tail()); //makes the first node of the original list point to the tail of the replacing list.
      (replacingseq->get_p_tail())->set_next(p_seek); //makes tail of replacing list join to original list
      original_list->set_p_head(replacingseq->get_p_head()); // makes the head of the joined lists the start of the replacing list.

      cout << "Next 10 base pairs: ";
      for (int i = 0; i < 10; i++) {
        if (p_seek != nullptr) {
          cout << p_seek->get_basepair();
          p_seek = p_seek->get_next();
        }
      }
    }

      cout << endl;



  }

  //if user wants to replace a part of the sequence that ends after the end of the sequence.
  else if ((bpposition + bplength) >= original_list->size()) {
    p_seek = original_list->get_p_tail();

    for (int i = 0; i < (original_list->size())-bpposition; i++) {
      p_seek = p_seek->get_previous(); //p_seek stops at the last node before the replacing sequence should begin.
    }
    p_temp = p_seek;
    cout << "Prev 10 base pairs: ";
    for (int i = 9; i >= 0; i--) {
      if (p_temp != nullptr) {
        previous[i] = p_temp->get_basepair();
        p_temp = p_temp->get_previous();
        num_previous++;
      }

    }
    for (int i = 10-num_previous; i < 10; i++) {
      cout << previous[i];
    }
    cout << endl;
    cout << "Region of interest: ";
    replacingseq->printlist();
    cout << endl;
    //connects the last node before replacement to the head of the replacing sequence.
    p_seek->set_next(replacingseq->get_p_head());
    (replacingseq->get_p_head())->set_previous(p_seek);
    original_list->set_p_tail(replacingseq->get_p_tail());
    cout << "Next 10 base pairs: ";
    if ((original_list->get_p_tail())->get_next() == nullptr) {
      cout << 0 << endl;
    } //There are no characters after the tail.

  }

  else {
    cout << "Prev 10 base pairs: ";
    for (int i = 0; i < bpposition-1; i++) {
      if ((bpposition - i) <= 10) {
        cout << p_seek->get_basepair();
      }
      p_seek = p_seek -> get_next(); //p_seek finishes at node before first node to be replaced.
    }
    cout << p_seek->get_basepair();
    cout << endl;
    cout << "Region of interest: ";
    replacingseq->printlist();
    cout << endl;
    for (int i = 0; i < bpposition + bplength; i++) {
      if (p_seek2 != nullptr) {
        p_seek2 = p_seek2 -> get_next(); // p_seek2 finishes at the node just after which the replacement should occur.
      }
    }
    p_temp = p_seek2; //to print next 10 pairs without changing p_seek2.
    cout << "Next 10 base pairs: ";
    for (int i = 0; i < 10; i++) {
      if (p_temp != nullptr) {
        cout << p_temp->get_basepair();
        p_temp = p_temp->get_next();
      }
    }

    cout << endl;
    //links the replacing seq with the original list after deleting the user specified length to be deleted in the original list.
    p_seek->set_next(replacingseq->get_p_head());
    (replacingseq->get_p_head())->set_previous(p_seek);
    p_seek2->set_previous(replacingseq->get_p_tail());
    (replacingseq->get_p_tail())->set_next(p_seek2);


 }

}

void Doublylinkedlist::search(Doublylinkedlist* desired_sequence, Doublylinkedlist* original_list) {
  Node* p_reference = original_list->get_p_head(); //keeps track of current position in the list.
  Node* p_matching = desired_sequence->get_p_head(); //keeps track of user searched list.
  Node* p_preceding; //for iterating to find previous 10 base pairs, likewise for p_next10.
  Node* p_next10;
  Node* p_temp; //temporary node pointer to take on the value of p_ref when the assessing whether there is a match.
  Node* p_iterator; //for printing out the region of interest.
  char previous[10];
  int last_index = 0, match_num = -1, count_previous = 0, count_loop = 0; //markers to count respective numbers. match_num is -1 since first match is match #0
  int original_size = original_list->size(), desired_size = desired_sequence->size(); //desired size is size of search pattern.

  //loops until the there can be no more matches.
  for (int i = 0; i <= (original_size - desired_size); i++) {
      last_index=i;
      p_temp = p_reference; //p_temp tries to see if there is a perfect match without changing p_reference (in case there is no match).
      p_iterator = p_reference;
      //this loop occurs to test whether there is a perfect match.
      while ((p_temp != nullptr) && (p_matching != nullptr) && (p_temp->get_basepair() == p_matching->get_basepair())) {

        count_loop++; //we only want p_preceding to be equal to the first matched base pair.
        if (count_loop == 1) {
        p_preceding = p_temp->get_previous(); //this node points to the first matched node
        }


        //condition of a complete match
        if (p_matching->get_next() == nullptr) {
          match_num++;
          last_index += desired_size-1;
          p_next10 = p_temp->get_next(); //p_next10 stores the location of the first node after

          cout << endl;
          cout << "Match Number #" << match_num << endl;
          cout << "Base pair positions: [" << i << ":" << last_index << "]" << endl;
          cout << "Base pair length: " << desired_size << endl;
          //if match is the first base pair.
          if (i == 0) {
            cout << "No preceding sequences available";
          }

          else {
            cout << "Prev 10 base pairs: ";
            for (int j = 9; j >= 0; j--) {
              if (p_preceding != nullptr) { //storing and printing previous 10 base pairs.
                previous[j] = p_preceding->get_basepair();
                p_preceding = p_preceding->get_previous();
                count_previous++;
              }

            }
            for (int j = 10-count_previous; j < 10; j++) {
              cout << previous[j];
            }
            count_previous = 0; //resets count previous to 0 for the next perfect match.
            }
            cout << endl;
            cout << "Region of interest: ";

            for (int k = 0; k < desired_size; k++) {
              cout << p_iterator->get_basepair();
              p_iterator = p_iterator->get_next();
            }
            cout << endl;
            //if match is at the end of the sequence.
            if (p_next10 == nullptr) {
              cout << "No pairs follow this sequence";
            }

            else {
              cout << "Next 10 base pairs: ";
              for (int j = 0; j < 10; j++) {
                if (p_next10 != nullptr) {
                  cout << p_next10->get_basepair();
                  p_next10 = p_next10->get_next();
                }
              }

            }
            cout << endl;


      }

      p_temp = p_temp->get_next();
      p_matching = p_matching->get_next();

    }
    //as long as after the match we have not gone past the last node.
    if (p_reference != nullptr) {
      count_loop = 0; //resets count loop.
      p_matching = desired_sequence->get_p_head(); //resets p_match to start from the first base pair.
      p_reference = p_reference -> get_next();

    }

  }

  if (match_num == -1) {
    cout << "No matches were found in this sequence." << endl;
  }

}

Doublylinkedlist::~Doublylinkedlist() {

  while (p_tail != nullptr) {
    pop_back();
  }
  cout << "deconstructor activated" << endl;
}

//stack that holds the nucleotide sequences from each file. (linked list of doublylinkedlist)
class Llofll {
private:
  Doublylinkedlist* p_headlist;
  Doublylinkedlist* p_taillist;
public:
  Llofll();
  void push(Doublylinkedlist* list_address);
  void pop();
  Doublylinkedlist* get_head();
  void print_listofl();
  ~Llofll();


};

Llofll::Llofll() {
  p_headlist = nullptr;
  p_taillist = nullptr;
}

void Llofll::push(Doublylinkedlist* newlist_address) {
  //if the stack is currently empty
  if (p_headlist == nullptr) {

    p_headlist = newlist_address;
    p_taillist = newlist_address;

  }

  else {

    Doublylinkedlist* p_temp = p_taillist;
    p_taillist->set_nextll(newlist_address);
    p_taillist = newlist_address;
    p_taillist->set_previousll(p_temp);

  }

}

void Llofll::pop() {
  Doublylinkedlist* p_temp;

  //if there is nothing inside the list
  if (p_taillist == nullptr) {
    cout << "there is nothing currently in the list to pop" << endl;
  }
  //if there is only one node remaining.
  else if (p_headlist == p_taillist) {
    p_temp = p_headlist;
    p_headlist = nullptr;
    p_taillist = nullptr;
    delete p_temp;
  }
  else {
    p_temp = p_taillist;
    p_taillist = p_taillist->get_previousll();
    p_taillist->set_nextll(nullptr);
    delete p_temp;
}
}

Doublylinkedlist* Llofll::get_head() {
  return p_headlist;
}
//prints everything in database.
void Llofll::print_listofl() {
  //goes to each list.
  for (Doublylinkedlist* p_seek = p_headlist; p_seek != nullptr; p_seek = p_seek->get_nextll()) {
    //goes to each node.
    for (Node* node_seek = p_seek->get_p_head(); node_seek != nullptr; node_seek = node_seek -> get_next()) {
      cout << node_seek->get_basepair();
    }
  }

}

Llofll::~Llofll() {
  while(p_taillist != nullptr) {
    pop();
  }
}

class DNADatabase {
private:
  int num_files = 2; //to store number of .fna files to open.
  Llofll sequences; //linked list of sequences
  vector<string> fna_files; //stores all the relevent file names, index corrresponds to sequence number in sequences linked list,
  int list_selection; //stores the index of the fna file we are interestd in.
  Doublylinkedlist* desired_sequence; //stores the user specified sequence
  int desired_position; //stores the user specified location for operation.
  Doublylinkedlist* sequence_tobeedited; //stores the latest version of the edited list.
  int bplength; //stores length of what the user wants for each function.

public:
  DNADatabase();
  bool check_ifinlist(string file); //checks if a specific file name is already stored in the list of file names that have already been loaded.
  void load_fna(vector <string>& files); //loads data from specified files into sequences member object and also stores appropriate file names in fna_files..
  bool check_if_fna (string filename); //returns 1 if file is an .fna file
  int display_userchoice(); //displays the files that are currently present in the data base for the user to make a selection.
  void set_list_selection(int choice); //takes in the users choice of which file to edit and stores it in the member variable choice
  void create_addedlist(char* user_pattern, int size); //transforms the user input to a linkedlist
  void get_list(); //isolates the list the user whats to edit by storing it in sequence_tobeedited.
  void insert_sequence(); //inserts the chosen sequence at a user defined point.
  void set_desiredposition(int position); //stores user specified position.
  void input_fromfile(string filename); //takes in input from an .fna file for options 2 and 3.
  void set_bplength(int length); //takes in user defined length and saves it to bplength.
  void delete_sequence(); //operations from options 2 and 3.
  void replace_sequence();
  void search_sequence();
  void save_tofile(string file);
  void search_database();
  ~DNADatabase();


};

DNADatabase::DNADatabase() {
  desired_sequence = nullptr;
}

bool DNADatabase::check_ifinlist(string file) {

  for (int i = 0; i < fna_files.size(); i++) {
    if (file == fna_files[i]) {
      return true;
    }
  }
  return false;
}

void DNADatabase::load_fna(vector <string>& files) {
  int size = files.size();
  int count = 1;
  for (int i = 0; i < size; i++) {

    if (check_if_fna(files[i]) && !check_ifinlist(files[i])) {
      //create fstream object and open the file
      fstream fnafile(files[i], fstream::in);
      if(!fnafile) {
        cout << "Failed to open " << files[i] << endl;
      }
      else {
        Doublylinkedlist* list = new Doublylinkedlist();
        fna_files.push_back(files[i]);
        cout << "Loading file #" << count << " \"" <<  files[i] << "\" ..." << endl;
        char c;
        if ((c = fnafile.get()) == '>') { //checks if there is a descriptor.
          do {
            c = fnafile.get(); //ignore the first line

          } while (c != '\n');
          do {
            if (c != '\n' && c != EOF) {
            list->append(c);

            }
            c = fnafile.get();
          } while (c != EOF);
        }
        else {

          do {
            if (c != '\n' && c != EOF) {
            list->append(c); //appends this the newly created list.

            }
            c = fnafile.get();
          } while (c != EOF);
          //ignore all newline characters.

        }

        sequences.push(list); //pushes list onto sequences which stores all the data.
        fnafile.close();
        count++;
      }

    }

  }
  cout << endl;
  for (int i = 0; i < size; i++) {
    files.pop_back(); //remove all file names so that files that were not approiate do not remain inside the vector should the user want to add more files.
  }
}

bool DNADatabase::check_if_fna (string filename) {
  if (filename[filename.size()-1] == 'a' && filename[filename.size()-2] == 'n' && filename[filename.size()-3] == 'f' && filename[filename.size()-4] == '.' ) {
    return true;
  }
  cout << filename << " is not a valid .fna file." << endl;
  return false;
}

int DNADatabase::display_userchoice() {
  cout << "(0) Return to main menu. " << endl;
  for (int i = 0; i < fna_files.size(); i++) {
    cout << "(" << i+1 << ") " << fna_files[i] << endl;
  }
  return (fna_files.size());
}

void DNADatabase::set_list_selection(int choice) {
  list_selection = choice;
}

void DNADatabase::create_addedlist(char* user_pattern, int size) {

  Doublylinkedlist* pattern = new Doublylinkedlist();
  for (int i = 0; i < size; i++) {
    pattern->append(user_pattern[i]);
  }
  desired_sequence = pattern; //desired sequence is the location of the linked list we want.

}

//must be called before an operation. Identifies the selected list based on users selection.
void DNADatabase::get_list() {
  Doublylinkedlist* p_seek = sequences.get_head();
  for (int i = 0; i < list_selection-1; i++) { //uses list_selection from private variables.
    p_seek = p_seek->get_nextll();
  }
  //p_seek is now the address of the linked list in database we are interested in.
  sequence_tobeedited = p_seek; //we set sequence_to be edited to be the list we want to edit.
}

void DNADatabase::insert_sequence() {
  while (desired_position > sequence_tobeedited->size()) {
    cout << "The position you entered is too long for the list." << endl;
    cout << "Please give a number that is smaller than " << sequence_tobeedited->size() << endl;
    cout << "> ";
    cin >> desired_position;
  }
  sequence_tobeedited->insert(sequence_tobeedited, desired_sequence, desired_position);

}

void DNADatabase::set_desiredposition(int position) {
  desired_position = position;
}

void DNADatabase::input_fromfile(string filename) {

  Doublylinkedlist* pattern = new Doublylinkedlist();
  desired_sequence = pattern;

  fstream fnafile(filename, fstream::in);
  if(!fnafile) {
    cout << "Failed to open " << filename << endl;
  }
  else {
    char c;
    if ((c = fnafile.get()) == '>') { //checks if there is a descriptor.
      do {
        c = fnafile.get(); //ignore the first line

      } while (c != '\n');
      do {
        if (c != '\n' && c != EOF) {
          desired_sequence->append(c);

        }
        c = fnafile.get();
      } while (c != EOF);
    }
    else {

      do {
        if (c != '\n' && c != EOF) {
          desired_sequence->append(c);

        }
        c = fnafile.get();

      } while (c != EOF);

    }

    fnafile.close();
  }

}

void DNADatabase::set_bplength(int length) {
  bplength = length;
}

void DNADatabase::delete_sequence() {
  while (desired_position >= sequence_tobeedited->size()) {
    cout << "The position you entered is too long for the list." << endl;
    cout << "Please give a number that is smaller or equal to " << sequence_tobeedited->size() << endl;
    cout << "> ";
    cin >> desired_position;
  }
  sequence_tobeedited->delete_seq(desired_position, bplength, sequence_tobeedited);
  }

void DNADatabase::replace_sequence() {
  while (desired_position >= sequence_tobeedited->size()) {
    cout << "The position you entered is too long for the list." << endl;
    cout << "Please give a number that is smaller than " << sequence_tobeedited->size() << endl;
    cout << ">";
    cin >> desired_position;
  }
  sequence_tobeedited->replace_seq(desired_position, bplength, sequence_tobeedited, desired_sequence);

}

void DNADatabase::search_sequence() {
  cout << "in dnadatabase function" << endl;
  sequence_tobeedited->search(desired_sequence, sequence_tobeedited);
}

void DNADatabase::save_tofile(string file) {
  ofstream newfile;
  newfile.open(file);
  for (Node* p_seek = sequence_tobeedited->get_p_head(); p_seek != nullptr; p_seek = p_seek->get_next()) {
    newfile << p_seek->get_basepair();
  }

  newfile.close();
  cout << endl;
  cout << "The file " << file << " has been successfully saved" << endl;
  cout << endl;
}

void DNADatabase::search_database() {
  int i = 0;
  for (Doublylinkedlist* p_seek = sequences.get_head(); p_seek != nullptr; p_seek = p_seek->get_nextll()) {
    //goes to each node.
    cout << "Analysing " << fna_files[i] << "..." << endl;
    p_seek->search(desired_sequence, p_seek);
    cout << endl;
    i++;

  }

}

DNADatabase::~DNADatabase() {
}

int main() {

  int user_choice1 = 0, user_choice2 = 0, user_choice3 = 0, size_filename, num_fnafilesinlist = 0, position = 0, length = 0;
  char current_char;
  vector <string> files_added; //stores all of the files given by user regardless of whether .fna.
  string to_beadded, desired_filename;
  char *sequence; //stores user inputted sequence.

  TicToc time;
  DNADatabase dna_db;

  cout << endl;
  cout << "Welcome to the DNA Editing program" << endl;


  while (user_choice1 != 4) {

    cout << endl;
    cout << "Select an option:" << endl;
    cout << "(1) Load DNA(s)." << endl;
    cout << "(2) Process a DNA." <<  endl;
    cout << "(3) Analyse the DNA database." << endl;
    cout << "(4) Quit." << endl;
    cout << ">";
    //if user puts an input that is not in the scope of 1 to 4.
    while (!(cin >> user_choice1) || user_choice1 > 4 || user_choice1 < 1) {
      cin.clear();
      while(cin.get() != '\n')
      continue;
      cout << "Please enter a number between 1 to 4" << endl;
      cout << ">";
    }

    cin.ignore(); // clears the newline character from the cin buffer.



    if (user_choice1 == 1) {
      cout << endl;
      cout << "Load DNA strands" << endl;
      cout << endl;
      cout << "Enter the DNA file names:" << endl;
      cout << "For multiple files, separate them by a comma. Only .fna are recognised." << endl;
      cout << ">";
      check_whitespace(files_added);
      cout << endl;
      time.tic();
      dna_db.load_fna(files_added);
      time.toc();
      cout << endl;
      cout << time << endl;

    }

    else if (user_choice1 == 2) {
      cout << endl;
      cout << "Select a DNA to process:" << endl;
      num_fnafilesinlist = dna_db.display_userchoice();

      while (!(cin >> user_choice2) || user_choice2 > num_fnafilesinlist || user_choice2 < 0) {
        cin.clear();
        while(cin.get() != '\n')
        continue;
        cout << "> Please enter a number from the selection" << endl;
      }

      if (user_choice2 != 0) {

        dna_db.set_list_selection(user_choice2);
        dna_db.get_list(); //identifies the selected list in our database.
        user_choice3 = 0;
        while (user_choice3 != 9) {

          cout << endl;
          cout << "Select from one of the following options" << endl;
          cout << "(1) Find DNA Sequence by input" << endl;
          cout << "(2) Find DNA sequence by file" << endl;
          cout << "(3) Add DNA sequence by input" << endl;
          cout << "(4) Add DNA sequence by file" << endl;
          cout << "(5) Delete DNA sequence by input" << endl;
          cout << "(6) Replace DNA sequence by input" << endl;
          cout << "(7) Replace DNA seqnece by file" << endl;
          cout << "(8) Save edited DNA sequence" << endl;
          cout << "(9) Exit submenu" << endl;
          cout << "> ";

          while (!(cin >> user_choice3) || user_choice3 < 0 || user_choice3 > 10) {
            cin.clear();
            while(cin.get() != '\n')
            continue;
            cout << "> Please enter a number from the selection" << endl;
          }


          if (user_choice3 == 1) {
            cout << endl;
            cout << "Find DNA sequence by input" << endl;
            cout << endl;
            cout << "Enter the DNA sequence to search (eg, GTCACT): " << endl;
            cout << "> ";

            cin >> to_beadded;
            time.tic();
            sequence = new char[to_beadded.size()];
            // converting string of input pattern to a dynamic character array.
            for (int i = 0; i < to_beadded.size(); i++) {
              sequence [i] = to_beadded[i];
            }
            dna_db.create_addedlist(sequence, to_beadded.size());
            dna_db.search_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;
            delete [] sequence;

          }
          else if (user_choice3 == 2) {
            cout << endl;
            cout << "Find DNA sequence by file" << endl;
            cout << endl;
            cout << "Enter the name of the .fna file containing the sequence to search:" << endl;
            cout << "> ";
            cin >> to_beadded;
            while (!dna_db.check_if_fna(to_beadded)) {
              cout << "The file you have entered is not of .fna extension. Please enter a valid file. " << endl;
              cout << "> ";
              cin >> to_beadded;
            }
            dna_db.input_fromfile(to_beadded);
            dna_db.set_desiredposition(position);
            time.tic();
            dna_db.search_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;
          }

          //add sequence by input
          else if (user_choice3 == 3) {
            cout << endl;
            cout << "Add DNA sequence by input" << endl;
            cout << endl;
            cout << "Enter a DNA sequence to add: " << endl;
            cout << "> ";
            cin >> to_beadded;
            sequence = new char[to_beadded.size()];
            // converting string of input pattern to a dynamic character array.
            for (int i = 0; i < to_beadded.size(); i++) {
              sequence [i] = to_beadded[i];
            }
            dna_db.create_addedlist(sequence, to_beadded.size());

            cout << "Enter a base pair position: " << endl;
            cout << "> ";
            cin >> position;
            dna_db.set_desiredposition(position);
            time.tic();
            dna_db.insert_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;
            delete [] sequence;

          }

          else if (user_choice3 == 4) {
            cout << endl;
            cout <<  "Add DNA sequence by file" << endl;
            cout << endl;
            cout << "Enter the name of the .fna file containing the sequence to add:" << endl;
            cout << "> ";
            cin >> to_beadded;
            while (!dna_db.check_if_fna(to_beadded)) {
              cout << "The file you have entered is not of .fna extension. Please enter a valid file. " << endl;
              cout << "> ";
              cin >> to_beadded;
            }
            dna_db.input_fromfile(to_beadded);
            cout << "Enter base pair position: " << endl;
            cout << "> ";
            cin >> position;
            dna_db.set_desiredposition(position);
            time.tic();
            dna_db.insert_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;

          }
          else if (user_choice3 == 5) {
            cout << endl;
            cout << "Delete DNA sequence by input" << endl;
            cout << endl;
            cout << "Enter a base pair position:" << endl;
            cout << "> ";
            cin >> position;
            dna_db.set_desiredposition(position);
            cout << "Enter a base pair length: " << endl;
            cout << "> ";
            cin >> length;
            dna_db.set_bplength(length);
            time.tic();
            dna_db.delete_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;


          }
          else if (user_choice3 == 6) {
            cout << endl;
            cout << "Replace DNA sequence by input" << endl;
            cout << endl;
            cout << "Enter nucleotide sequence you would like to add: " << endl;
            cout << "> ";
            cin >> to_beadded;

            sequence = new char[to_beadded.size()];
            // converting string of input pattern to a dynamic character array.
            for (int i = 0; i < to_beadded.size(); i++) {
              sequence [i] = to_beadded[i];
            }
            dna_db.create_addedlist(sequence, to_beadded.size());
            cout << "Enter a base pair position: " << endl;
            cout << "> ";
            cin >> position;

            dna_db.set_desiredposition(position);

            cout << "Enter a base pair length" << endl;
            cout << "> ";
            cin >> length;

            dna_db.set_bplength(length);
            time.tic();
            dna_db.replace_sequence();
            time.toc();
            cout << endl;
            cout << time << endl;

            delete [] sequence;
          }
          else if (user_choice3 == 7) {
            cout <<  endl;
            cout << "Replace DNA sequence by file" << endl;
            cout << endl;

            cout << "Enter the name of the .fna file containing the sequence you wish to replace into the file:" << endl;
            cout << "> ";
            cin >> to_beadded;
            while (!dna_db.check_if_fna(to_beadded)) {
              cout << "The file you have entered is not of .fna extension. Please enter a valid file. " << endl;
              cout << "> ";
              cin >> to_beadded;
            }
            dna_db.input_fromfile(to_beadded);
            cout << endl;
            cout << "Enter a base pair position: " << endl;
            cout << "> ";
            cin >> position;

            dna_db.set_desiredposition(position);
            cout << endl;
            cout << "Enter base pair length to be replaced" << endl;
            cout << "> ";
            cin >> length;

            dna_db.set_bplength(length);

            time.tic();
            dna_db.replace_sequence();
            time.toc();
          }
          else if (user_choice3 == 8) {
            cout <<  endl;
            cout << "Save edited DNA sequence" << endl;
            cout << endl;
            cout << "Enter a name you would like the file to be saved as (with .fna extension)" << endl;
            cout << "> ";
            cin.ignore();
            getline(cin, desired_filename);
            time.tic();
            dna_db.save_tofile(desired_filename);
            time.toc();
            cout << endl;
            cout << time << endl;
          }

        }
      }

    }

    else if (user_choice1 == 3) {
      cout << endl;
      cout << "Analyse DNA Database" << endl;
      cout << endl;
      cout << "Enter the name of the .fna file containing the sequence to search:" << endl;
      cout << "> ";
      cin >> to_beadded;
      while (!dna_db.check_if_fna(to_beadded)) {
        cout << "The file you have entered is not of .fna extension. Please enter a valid file. " << endl;
        cout << "> ";
        cin >> to_beadded;
      }
      dna_db.input_fromfile(to_beadded);
      time.tic();
      dna_db.search_database();
      time.toc();
      cout << endl;
      cout << time << endl;

    }

  }

  cout << endl;
  cout << "Exiting program..." << endl;
return 0;
}
//checks user input for file names for white spaces and removes leading and ending spaces.
void check_whitespace(vector<string> &filenames) {
  int char_flag = 0, space_flag = 0, sizeoffile_name = 0; //to keep track of whether the white space we are at is necessary.
  char current_char; //stores value of character we are currently looking at.
  string file_name; //stores the correct file name after whitespace analysis.
  do {
    current_char = cin.get();
    if (current_char != ' ' && current_char != ',' && current_char != '\n') {
      char_flag++;
      file_name.push_back(current_char);
      space_flag = 0; //resets space_flag.
    }
    else if (char_flag > 0 && current_char == ' ') {
      space_flag++;
      file_name.push_back(current_char);
    }
    else if (current_char == ',' || current_char == '\n') {
      char_flag = 0; //resets char flag.

      for (int i = 0; i < space_flag; i++) {
        file_name.pop_back(); //removes the end white spaces.
      }
      space_flag = 0;
      filenames.push_back(file_name);
      sizeoffile_name = file_name.size();
      for (int i = 0; i < sizeoffile_name; i++) {

        file_name.pop_back(); //removes all characters inside file_name so that a new file_name can generated.
      }

    }

  } while(current_char != '\n');

}
