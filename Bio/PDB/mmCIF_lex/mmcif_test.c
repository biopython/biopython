#include <stdio.h>

main()
{
  int flag=1;

  mmcif_set_file(stdin);

  while(flag)
    {
      flag=mmcif_get_token();
      if(!flag)
	{
	  return;
	}	
      printf("%d ", flag);
      printf("%s\n", mmcif_get_string());
    }	
}	

