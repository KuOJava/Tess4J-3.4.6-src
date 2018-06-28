import java.awt.Component;
import java.awt.Container;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Label;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.IOException;

import javax.security.auth.login.CredentialExpiredException;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

//������������ͼ���û������д������
public class Main {
	private JTextField inputPath;
	private JButton choosePath;
	private JButton confirm;
	private JLabel situation;
	public static JTextField showMessage;
	private JLabel jumpToExcel;
	private JButton clickToCopy;

	public Main() {
		// ����һ�������Ϊ��������
		// JFrame.setDefaultLookAndFeelDecorated(true);
		JFrame f = new JFrame("����ͼƬʶ��ϵͳ");
		// f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		// ��ȡ����е��������
		Container contentPane = f.getContentPane();
		// �������񲼾�
		contentPane.setLayout(new GridBagLayout());
		// ʹ�����񲼾�Լ����
		GridBagConstraints c = new GridBagConstraints();
		// �������,���������ĵ����ṩ��ͼ�ν����
		// '����·��'�ı��򣺴��ڵ�1�У�ռ3��;'ѡ��·��'��ť�����ڵ�1�У���4�У�'ȷ��'��ť�����ڵ�1�У���5��
		// 'ʶ�����'��ʾ�򣺴��ڵ�2�У�ռ������
		// ��ʾ��Ϣ���ı��򣺴��ڵ�3,4�У�ռ������
		// '������������Excel�ĸ�Ŀ¼'�����ӣ����ڵ�5�У�ռ3�У�'��������ļ�'��ť�����ڵ�5�У�ռ����

		// '����·��'���ı���
		inputPath = new JTextField(1);
		inputPath.setText("././b");
		inputPath.setFont(new Font("����", Font.PLAIN, 80));
		inputPath.setHorizontalAlignment(JTextField.CENTER);
		// ��һ�е�һ�п�ʼ
		c.gridx = 0;
		c.gridy = 0;
		// ��3��
		c.gridwidth = 3;
		c.gridheight = 1;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 3;
		c.weighty = 1;
		// ���õ���¼�

		// ����c��������������
		contentPane.add(inputPath, c);

		// 'ѡ��·��'��ť
		choosePath = new JButton("ѡ��·��");
		choosePath.setFont(new Font("����", Font.PLAIN, 50));
		// ѡ��·������¼�
		choosePath.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				Component parent = null;
				int returnVal = chooser.showOpenDialog(parent);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					inputPath.setText(chooser.getSelectedFile().getAbsolutePath());
				}
			}
		});
		// ���ڵ�һ�е�����
		c.gridx = 3;
		c.gridy = 0;
		c.gridwidth = 1;
		c.weightx = 1;
		c.weighty = 1;
		// c.gridheight=1;----���ø�ֵ
		contentPane.add(choosePath, c);
		// 'ȷ��'��ť
		confirm = new JButton("ȷ��");
		confirm.setFont(new Font("����", Font.PLAIN, 50));
		// ȷ����ť�ĵ���¼�
		confirm.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// ������ȡ·����ȡ����
				Long start = System.nanoTime();
				String path = inputPath.getText();
				if (path.equals("")) {
					showMessage.setText("��û�������κ�·��\n");
				} else {
					showMessage.setText("����ʶ�𡣡���������");
					Thread thread =new Thread(new Runnable() {
						@Override
						public void run() {
							String message = "";
							try {
								message = Tess.find(path, Tess.depth);
							} catch (IOException e) {
								e.printStackTrace();
							}
							PoiExcel poiexcel = new PoiExcel();
							poiexcel.poi(Tess.temp1, Tess.temp2);// ����Excel���
							showMessage.setText(message);
							if(message.equals("���")) {
								Long end = System.nanoTime();
								System.out.println("ʱ�䣺"+(end-start)+"����");
								System.out.println("�ٶȣ�"+(end-start)/Tess.count+"����/��"+"\t��\t"+(end-start)/Tess.count/1000000+"����/��");
							}
						}
					});
					thread.start();
				}
			}
		});
		c.gridx = 4;
		// c.gridy=4;----���ø�ֵ
		// c.gridwidth=1;----���ø�ֵ
		// c.gridheight=1;----���ø�ֵ
		c.weightx = 1;
		c.weighty = 1;
		contentPane.add(confirm, c);

		// 'ʶ�����:'��ʾ��
		situation = new JLabel("ʶ�����:");
		situation.setFont(new Font("����", Font.PLAIN, 50));
		situation.setHorizontalAlignment(SwingConstants.CENTER);
		c.gridx = 0;
		c.gridy = 1;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth=5;
		c.gridheight=1;
		contentPane.add(situation, c);

		// '��ʾ��Ϣ�Ŀ�'
		showMessage = new JTextField("");
		showMessage.setFont(new Font("����", Font.PLAIN, 50));
		showMessage.setHorizontalAlignment(JTextField.CENTER);
		// c.gridx=0;
		c.gridy = 2;
		c.gridwidth = 5;
		c.gridheight = 2;
		c.weightx = 5;
		c.weighty = 2;
		contentPane.add(showMessage, c);

		// '������������Excel�ĸ�Ŀ¼'
		jumpToExcel = new JLabel("����������Excel:");
		jumpToExcel.setFont(new Font("����", Font.PLAIN, 50));
		jumpToExcel.setHorizontalAlignment(SwingConstants.CENTER);
		// ����¼�
		c.gridx = 0;
		c.gridy = 4;
		c.gridwidth = 3;
		c.gridheight = 1;
		c.weightx = 3;
		c.weighty = 1;
		contentPane.add(jumpToExcel, c);

		// '��������ļ�'
		clickToCopy = new JButton("���ļ�");
		clickToCopy.setFont(new Font("����", Font.PLAIN, 50));
		// ���ļ�
		clickToCopy.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				StartExcel startexcel = new StartExcel();
				startexcel.start();// �Զ���Excel���
			}
		});
		c.gridx = 3;
		c.gridwidth = 2;
		c.weightx = 2;
		c.weighty = 1;
		contentPane.add(clickToCopy, c);

		// ��ʾ����
		f.pack();
		//������ʾ
		int width =1200;
		int height = 900;
		Point point = GraphicsEnvironment.getLocalGraphicsEnvironment().getCenterPoint();
		f.setBounds(point.x-width/2, point.y-height/2, width, height);
		//f.setSize(1200, 900);
		f.setVisible(true);
		// ���ô��ڵĹر�
		f.addWindowListener(new WindowListener() {

			@Override
			public void windowOpened(WindowEvent e) {
				// TODO Auto-generated method stub

			}

			@Override
			public void windowIconified(WindowEvent e) {
				// TODO Auto-generated method stub

			}

			@Override
			public void windowDeiconified(WindowEvent e) {
				// TODO Auto-generated method stub

			}

			@Override
			public void windowDeactivated(WindowEvent e) {
				// TODO Auto-generated method stub

			}

			@Override
			public void windowClosing(WindowEvent e) {
				// TODO Auto-generated method stub
				System.exit(0);
			}

			@Override
			public void windowClosed(WindowEvent e) {
				// TODO Auto-generated method stub

			}

			@Override
			public void windowActivated(WindowEvent e) {
				// TODO Auto-generated method stub

			}
		});
	}

	public static void main(String[] args) {
		new Main();
	}
}
