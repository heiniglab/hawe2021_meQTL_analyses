BEGIN { n=split(key, order, ":") }
{
   for(i=1; i<n; i++) {
     printf("%s\t", $order[i])
    }
    printf("%s", $order[n])
    printf("\n");
}
