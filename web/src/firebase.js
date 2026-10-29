import { initializeApp } from "firebase/app";
import { getAuth, GoogleAuthProvider } from "firebase/auth";

// PanoptesBio Web App Firebase configuration
export const firebaseConfig = {
  apiKey: "AIzaSyB4nBuxqnKFqFeyS_RZVCNlK6ox0iXEbq8",
  authDomain: "panoptesbio.firebaseapp.com",
  projectId: "panoptesbio",
  storageBucket: "panoptesbio.firebasestorage.app",
  messagingSenderId: "561444549088",
  appId: "1:561444549088:web:e1be29926bdaaf9766fb85",
  measurementId: "G-3CQM0RGVH3"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
export const auth = getAuth(app);
export const googleProvider = new GoogleAuthProvider();

export default app;
