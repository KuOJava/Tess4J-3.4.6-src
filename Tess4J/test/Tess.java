import java.io.File;  
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import LSD.LSD;
import LSD.Line;
public class Tess {  
	 static String [] retainWord = new String[]{"��ҵע���","��ҵ����","����","����","ǩ֤����","У׼ʱ��","��Ӫ��Χ","����ʱ��","ע���ʱ�","ס��"};//������ 
     static ArrayList<String> temp1=new   ArrayList<String>();
     static ArrayList<String> temp2=new   ArrayList<String>();
     static boolean tem1=true;
     static boolean tem2=true;
     static boolean tem3=false;
     static boolean tem4=false;
     static  String h=null;
  
     static int depth=1;  
     static int count = 0;
    
  /***����ָ���ļ���ȡ��ȡͼƬ****/    
    public static String find(String pathName,int depth) throws IOException{  
        //��ȡpathName��File����  
        File dirFile = new File(pathName);  
        //�жϸ��ļ���Ŀ¼�Ƿ���ڣ�������ʱ�ڿ���̨�������  
        if (!dirFile.exists()) {  
            System.out.println("do not exit");  
            return "�����ڸ��ļ���";  
        }  
        //�ж��������һ��Ŀ¼�����ж��ǲ���һ���ļ���ʱ�ļ�������ļ�·��  
        if (!dirFile.isDirectory()) {  
            if (dirFile.isFile()) {  
                System.out.println(dirFile.getCanonicalFile());  
            }  
            return "�����ļ���";  
        }  
        //��ȡ��Ŀ¼�µ������ļ�����Ŀ¼��  
        String[] fileList = dirFile.list();  
        int currentDepth=depth+1;  
        for (int i = 0; i < fileList.length; i++) {  
            //�����ļ�Ŀ¼  
            String string = fileList[i];  
            //File("documentName","fileName")��File����һ��������  
            File file = new File(dirFile.getPath(),string);  
            String name = file.getName();  
            //�����һ��Ŀ¼���������depth++�����Ŀ¼���󣬽��еݹ�  
            if (file.isDirectory()) {  
                //�ݹ�  
                find(file.getCanonicalPath(),currentDepth);  
            }else{  
                //�ж��Ƿ���ͼƬ��ʽ������ǵĻ�������ȡ��������ݣ���������ҵ���ƺ�ע��ŵĶ�Ӧ���ݴ��뵽��Ӧ���б���
            	if (name.endsWith(".png")||name.endsWith(".jpg")) 
            	{
            		File imageFile = new File(dirFile+"\\"+name);
//                   System.out.println(dirFile+"\\"+name);  
            		BufferedImage img = ImageIO.read(imageFile);
            		//������������Ƭ�Ļ�ִ�и�����--lsz
            		if(name.endsWith(".jpg")) {
            			double max_x1=0;
            			double max_x2=0;
            			double max_y1=0;
            			double max_y2=0;
            			
            			double min_x1=0;
            			double min_x2=0;
            			double min_y1=0;
            			double min_y2=0;

            			double min_x11=0;
            			double min_x22=0;
            			double min_y11=0;
            			double min_y22=0;
            			
            			double max_x11=0;
            			double max_x22=0;
            			double max_y11=0;
            			double max_y22=0;
            			Graphics2D g2d = img.createGraphics();
            			int x = img.getWidth();
            			int y = img.getHeight();
            			max_x1=x;
            			//System.out.println(x);
            			//System.out.println(y);
            		
            			HashSet<Line> lines = new HashSet<Line>();

            			
            			double [] arr = img.getData().getPixels(0,0,x,y,new double[x*y*3]);

            			double [] arr2 = new double[x*y];
            		
            			//System.out.println("arr.length:"+arr.length);
            			int c=0;
            			for(int i1 = 0; i1 < arr.length-3; i1+=3) {
            				double B = arr[i1];
            				double G = arr[i1+1];
            				double R = arr[i1+2];
            				double level = R * 0.2126 + G * 0.7152 + B * 0.0722;
            				arr2[c++] = level;
            			}

            			LSD lsd = new LSD();

            			double [] out = lsd.lsd(arr2,x,y);

            			for(int i1 = 0; i1 < lsd.n_out; i1++) {
            				for (int j = 0; j < 7; j++)
            				
            				lines.add(new Line(out[7 * i1 + 0], out[7 * i1 + 1],
            						out[7 * i1 + 2], out[7 * i1 + 3]));

            			}
            		
            			
            			int z=0;
            			//System.out.println("��ʼ");
            			for (Line l : lines) {
            				if((l.x2-l.x1>x/8)&&(l.y1>min_y11)&&(l.y1<y/2)&&(Math.abs(l.y1-l.y2)<10))
            				{
            					//System.out.println("true");
            					min_x11=l.x1;
            					min_x22=l.x2;
            					min_y11=l.y1;
            					min_y22=l.y2;
            				}
            				if((l.x2-l.x1>x/8)&&(l.y1<max_y11)&&(l.y1>y/2)&&(Math.abs(l.y1-l.y2)<10))
            				{
            					//System.out.println("true");
            					max_x11=l.x1;
            					max_x22=l.x2;
            					max_y11=l.y1;
            					max_y22=l.y2;
            				}
            				
            				
            				if(l.x2-l.x1>x/8)
            				{	
            					++z;
//            				System.out.println(l.y2-l.y1);
//            				System.out.println(l.x1);
//            				System.out.println(l.x2);
//            				System.out.println(l.y1);
//            				System.out.println(l.y2);
//            				System.out.println(z);
            				}
            				if((l.y2-l.y1>y/8)&&(l.x1>min_x1)&&(l.x1<x/2)&&(Math.abs(l.x1-l.x2)<10))
            				{
            					min_x1=l.x1;
            					min_x2=l.x2;
            					min_y1=l.y1;
            					min_y2=l.y2;
            				}
            				else if((l.y2-l.y1>y/8)&&(l.x1<max_x1)&&(l.x1>x/2)&&(Math.abs(l.x1-l.x2)<10))
            				{
            					
            					max_x1=l.x1;
            					max_x2=l.x2;
            					max_y1=l.y1;
            					max_y2=l.y2;
            				}
            				
            			}
            			g2d.setPaint(Color.red);
            		//	g2d.drawLine((int)max_x1,(int)max_y1,(int)max_x2,(int)max_y2);
            		
            		
            			//g2d.drawLine((int)min_x1,(int)min_y1,(int)min_x2,(int)min_y2);
            			g2d.drawLine((int)min_x11,(int)min_y11,(int)min_x22,(int)min_y22);
            		
//            			System.out.println(max_x1+" "+min_x1+" "+min_y11);
//            			System.out.println(x/2);
//            			System.out.println(y/8);
//            			System.out.println("����");
//            			System.out.println(lines.size());
            			img=img.getSubimage((int)min_x1,(int)min_y11,(int)(max_x1-min_x1),(int)(y-min_y11));
            			//��ʾ���и��ͼƬ
            			JFrame jFrame = new JFrame("");
            			jFrame.add(new JLabel(new ImageIcon(img)));
            			jFrame.setVisible(true);
            		}
            		shadow shadows=new shadow();
            		try {    
            			String result =shadows.ocr(img) ;            			
            			String [] classify0=result.split("\n");
            			String [] classify=new String[classify0.length];
            			for(int j=0;j<classify0.length;j++) {
            				classify[j]="";
            				char[] a=classify0[j].toCharArray();
            				for(int k=0;k<classify0[j].length();k++) {
            					if(a[k]!=' ') {
            						classify[j]+=a[k];
            					}
            				}
            			}
            			for(int i1=0;i1<classify.length;i1++)
            			{
            				if(classify[i1].toString().indexOf("��ҵע���")!=-1||classify[i1].toString().indexOf("��ҵ����")!=-1)
            				{
            					if(classify[i1].toString().indexOf("��ҵ����")!=-1)
            					{   
            						if(classify[i1].toString().indexOf("��")!=-1) {
            							h=classify[i1].substring(classify[i1].indexOf("��")+1 );
                						temp1.add(h);
                						tem1=false;
            						}
            						else {
            							h=classify[i1].substring(classify[i1].indexOf("��")+1 );
                						temp1.add(h);
                						tem1=false;
            						}
            					}
            					else
            					{
            						if(classify[i1].toString().indexOf("��")!=-1) {
            							h=classify[i1].substring(classify[i1].indexOf("��")+1 );
            							temp2.add(h);
            							tem2=false;
            						}
            						else {
            							h=classify[i1].substring(classify[i1].indexOf("��")+1 );
            							temp2.add(h);
            							tem2=false;
            						}
            					}
            					
            				}
            				
            			}
            			System.out.println(result);
            			if(tem1)
            				temp1.add("none");
            			if(tem2)
            				temp2.add("none");
            			//	System.out.println(classify[i]);     
            			Main.showMessage.setText("�Ѿ�ʶ���:"+(++count)+"��");
            		} catch (Exception e) {    
            			//System.err.println(e.getMessage());    
            		} 
            	}
            }  
        }  
        return "���";
    }  
      
   /* public static void main(String[] args) throws IOException{  
        find("././b", depth);  //����ͼƬ����ȡָ�����ݣ�������ͼƬ���λ�õĴ洢·������Ҫ�޸�
  	    PoiExcel poiexcel=new PoiExcel();
	    poiexcel.poi(temp1,temp2);//����Excel���
		StartExcel startexcel=new StartExcel();
	    startexcel.start();//�Զ���Excel���
    }  */
}  