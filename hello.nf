#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = Channel.of(params.greeting)

process SPLITLETTERS {
    input:
    val x //the variable x will be taking inputs from our channel greeting_ch

    output:
    path 'chunk_*'
    //the output is a set of files with the names chunk_*
    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch) //the output of SPLITLETTER is a set of files, that we turn into a channel.
    results_ch = CONVERTTOUPPER(letters_ch.flatten()) // The .flatten transforms every element of the letters channel into its own item.
    results_ch.view{ it } //
}

