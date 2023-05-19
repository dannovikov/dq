process first_process {
    output:
    stdout()

    script:
    """
    echo "Hello, World!"
    """
}

process second_process {
    input:
    stdin

    output:
    path("output.txt")

    script:
    """
    cat - > output.txt
    """
}

workflow {
    my_output = first_process()
    second_process(my_output)
}
