//                  -*- mode:c++; tab-width:4 -*-
// File: tools.cpp ----------------------------------------------------------
// The C++ tools library implementation files.
// Assorted utility routines for use in C++ programs.
// Authors:  Alice E. Fischer and Michael J. Fischer.
// Modified September 2012; September 2016.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#include "tools.hpp"

// --------------------------------------------------------------------------
// I/O Utilities. -----------------------------------------------------------
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Used to flush the rest of the input line as in cin >> x >> flush;
// or cin >> flush;
istream&
flush(istream& is) {
    ios::iostate save = is.rdstate();
    is.clear();
    is.ignore(numeric_limits<streamsize>::max(), '\n');
    is.setstate( save );
    return is;
}

// --------------------------------------------------------------------------
// Print a neat header on the output.

void fbanner(ostream& fout) {
        char date[16], time[10];
        when(date, time);
        fout
                        << "\n---------------------------------------------------------------\n"
                        << "\t" << NAME << "\n\t" << CLASS << "\n\t" << date << "\t" << time
                        << endl
                        << "---------------------------------------------------------------\n";
}

// ---------------------------------------------------------------------------
//  Print termination message.

void bye(void) {
        cout << "\n---------------------------------------------------------------"
                        "\nNormal termination.\n";
}

// ----------------------------------------------------------------------------
//  Error handling and error recovery.
//-----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// This function is for error messages.
//    It takes a format argument followed by any number of data arguments.
//    It formats and prints an error message, then exits.

void
fatal (const char* format, ...)
{   va_list vargs;                               // optional arguments

    cout << flush;
    va_start(vargs, format);
    vfprintf(stdout, format, vargs);
    fprintf(stdout, "\n");
    va_end(vargs);
    exit(1);
}

// ----------------------------------------------------------------------------
// This is the constructor for the Fatal exception class.
//    It takes a format argument followed by any number of data arguments.
//    It formats and stores an error message.
//    The stored string is what C's printf function would print given
//    the same arguments.

Fatal::Fatal(const char* format, ...) {
        va_list vargs;                               // optional arguments
        va_start(vargs, format);
        vsnprintf(msg, 256, format, vargs);
}

// ----------------------------------------------------------------------------
// Routines for handling the time and date.
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Store the current date and time in the arguments.
//      System's date format is: "Fri Jun  9 10:15:55 1995\n"
//      After extraction, date is: "Fri Jun  9 1995"    hour is: "10:15:55"

void when(char* date, char* hour) {
        time_t now;           // Stores an integer encoding of the date and time.
        char* nowstring;    // Stores the date and time in a readable form.

        now = time(nullptr);              // Get the date and time from the system.
        nowstring = ctime(&now);                   // Convert to string form.
        strncpy(date, nowstring, 10);             // Extract day, month, date.
        strncpy(&date[10], &nowstring[19], 5);    // Extract space and year.
        date[15] = '\0';                          // Add the string terminator.
        strncpy(hour, &nowstring[11], 8);         // Copy hour:minutes:seconds.
        hour[8] = '\0';                           // Add the string terminator.
}

// ----------------------------------------------------------------------------
// Store the current date in the argument and return a pointer to it.
//      date format is: "Fri Jun  9 1995"

char* today(char* date) {
        time_t now;          // Stores an integer encoding of the date and time.
        char* nowstring;   // Stores the date and time in a readable form.

        now = time(nullptr);              // Get the date and time from the system.
        nowstring = ctime(&now);                   // Convert to string form.
        strncpy(date, nowstring, 10);             // Extract day, month, date.
        strncpy(&date[10], &nowstring[19], 5);    // Extract space and year.
        date[15] = '\0';                          // Add the string terminator.
        return (date);
}

// ----------------------------------------------------------------------------
// Store the current time in hour and return a pointer to it.
//      hour format is: "10:15:55"

char* oclock(char* hour) {
        time_t now;          // Stores an integer encoding of the date and time.
        char* nowstring;   // Stores the date and time in a readable form.

        now = time(nullptr);              // Get the date and time from the system.
        nowstring = ctime(&now);              // Convert to string form.
        strncpy(hour, &nowstring[11], 8);    // Extract hour, minutes, seconds.
        hour[8] = '\0';                      // Add the string terminator.
        return (hour);
}
